using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.Direct;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using TriangleNet.Meshing.Algorithm;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled7and9eqsSolutionex6
    {
        const double Sc = 0.1;

        private const double timeStep = 0.001; // in sec
        const double totalTime = 0.04; // in sec
        static int incrementsPertimeStep = 1;


        // strucutral model Loads
        static StructuralDof loadedDof = StructuralDof.TranslationZ;
        static double load_value = 0; //[kN]


        //structural model properties
        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa

        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa

        static int currentTimeStep = 0;
        static double lambda0 = 1;
        int nGaussPoints = 1;
        static double initial_dp_dx = 0.0;
        static double initial_dp_dy = 0.0;
        static double initial_dp_dz = 0.0;
        static double velocityDivInitialVal = 0;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>();
        private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();


        static double pressureMonitorNodeX = 0.05;
        static double pressureMonitorNodeY = 0.05;
        static double pressureMonitorNodeZ = 0.05;

        static double[] pressureMonitorNode = new double[]
            { pressureMonitorNodeX, pressureMonitorNodeY, pressureMonitorNodeZ };

        private static int pressureMonitorID;

        static double divPMonitorGPX = 0.01714618146;
        static double divPMonitorGPY = 0.05114110504;
        static double divPMonitorGPZ = 0.06533804930;

        //static double divPMonitorGPX = 0.05;
        //static double divPMonitorGPY = 0.05;
        //static double divPMonitorGPZ = 0.05;

        static double[] divPMonitorGP = new double[] { divPMonitorGPX, divPMonitorGPY, divPMonitorGPZ };

        private static int divPMonitorID;

        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationX;

        //Darcy model properties

        static double k_th = 7.52e-10; // m2/(KPa sec)
        static double Sv = 7e+3; // 1/(m)
        static double Lp = 2.7e-12; // m/(KPa sec)
        static double LplSvl = 0; // 1/(KPa sec)
        static double pv = 4; // kPa
        static double pl = 0d; // KPa
        static double div_vs = 1e-6; // 1/(sec)

        //static int nodeIdToMonitor = 36;

        //The above coordinates are in the reference configuration.
        //The actual coordinates will be computed in the analysis.

        static double structuralMonitorNodeX = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        static double structuralMonitorNodeY = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        static double structuralMonitorNodeZ = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa

        static double[] structuralMonitorNode = new double[]
            { structuralMonitorNodeX, structuralMonitorNodeY, structuralMonitorNodeZ };

        private static int structuralMonitorID;
        static ConvectionDiffusionDof eq7n8dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;

        //private static double boundaryValueAllBoundaries = pv;
        //private static double initialCondition = pv;
        private static double boundaryValueAllBoundaries = 0;
        private static double initialCondition = 0;

        static double modelMinX = 0;
        static double modelMaxX = 0.1;
        static double modelMinY = 0;
        static double modelMaxY = 0.1;
        static double modelMinZ = 0;
        static double modelMaxZ = 0.1;

        public Coupled7and9eqsSolutionex6()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]
        //[InlineData("../../../DataFiles/workingTetMesh4886.mphtxt")]
        //[InlineData("../../../DataFiles/chipMelter2M.mphtxt")]
        //[InlineData("../../../DataFiles/MeshCyprusTM.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh4886_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh648_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh648_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingQuadMesh64_1Domain.mphtxt")]
        [InlineData("../../../DataFiles/workingQHexaMesh6x6x6_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingQuadMesh27_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh155.mphtxt")]
        public void MonophasicEquationModel(string fileName)
        {
            //ContinuumElement3DGrowth.dT = timeStep;

            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                lambda.Add(elem.Key, lambda0);
            }

            Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var gpTensorDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    gpTensorDiv[i1] = new double[] { initial_dp_dx, initial_dp_dy, initial_dp_dz };
                }

                pressureTensorDivergenceAtElementGaussPoints.Add(elem.Key, gpTensorDiv);
            }

            Dictionary<int, double[]> velocityDivergenceAtElementGaussPoints =
                new Dictionary<int, double[]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocityDiv[i1] = velocityDivInitialVal;
                }

                velocityDivergenceAtElementGaussPoints.Add(elem.Key, velocityDiv);
            }

            miNormal = miTumor;
            kappaNormal = kappaTumor;
            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNode, 1e-2);
            //pressureMonitorID = Utilities.FindRandomInternalNode(comsolReader.NodesDictionary, modelMinX, modelMaxX,
            //modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            pressureMonitorID = structuralMonitorID;
            var eq78Model = new Eq78ModelProviderForStaggeredSolutionex6(comsolReader, k_th, Lp, Sv, pv, LplSvl, pl,
                velocityDivergenceAtElementGaussPoints,

                boundaryValueAllBoundaries, initialCondition, pressureMonitorID, eq7n8dofTypeToMonitor, modelMinX,
                modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);

            var eq9Model = new Eq9ModelProviderForStaggeredSolutionex6(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, eq9dofTypeToMonitor, loadedDof, load_value, modelMinX, modelMaxX, modelMinY,
                modelMaxY, modelMinZ, modelMaxZ);

            var equationModel = new Coupled7and9eqsModelex6(eq78Model, eq9Model, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, timeStep,
                totalTime, incrementsPertimeStep);



            #region loggin

            var p_i = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);

            double[] modelMaxVelDivOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dxOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dyOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dzOverTime = new double[(int)(totalTime / timeStep)];
            //var dp_dxi = new List<double[]>();
            //dp_dxi.Add(modelMax_dP_dxOverTime);
            //dp_dxi.Add(modelMax_dP_dyOverTime);
            //dp_dxi.Add(modelMax_dP_dzOverTime);
            #endregion


            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging


                var divPMonitorId =
                    Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], divPMonitorGP, 1e-1);

                //nodal logs
                p_i[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), ConvectionDiffusionDof.UnknownVariable];
                double[] allValues = new double[] { p_i[currentTimeStep] };
                structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];



                //modelMax_dP_dxOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .xcoeff_OverTimeAtGp1[0];
                //modelMax_dP_dyOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .ycoeff_OverTimeAtGp1[0];
                //modelMax_dP_dzOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .zcoeff_OverTimeAtGp1[0];



                modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints
                    .Select(x => Math.Abs(x.Value[0])).ToArray().Max();

                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints
                //    .Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }

                #endregion

                (equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;

            //var outputPath = $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\";
            var outputPath = "../../../StaggeredSolutionPresDynamex6/results1/";
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsX, outputPath + $@"pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{1}_.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsY, outputPath + $@"pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{2}_.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsZ, outputPath + $@"pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{3}_.txt");
            //(new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(modelMaxVelDivOverTime, outputPath + $@"pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_modelMaxVelDivOverTime.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsZ, outputPath + $@"pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_modelMax_dP_dx{3}OverTime.txt");


            var path = outputPath + "dp_dxi_mslv.csv";
            //CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredSolutionPresDynamex6/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredSolutionPresDynamex6/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredSolutionPresDynamex6/pi_nodes_mslv.csv");

        }






    }
}
