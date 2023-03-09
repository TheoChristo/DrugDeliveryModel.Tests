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
using System.Xml.Linq;
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class UPTdummySolution
    {
        const double Sc = 0.1;

        private const double timeStep = 0.00001; // in sec
        const double totalTime = 0.0001; // in sec
        static int incrementsPertimeStep = 1;
        static int currentTimeStep = 0;

        #region Structural model properties

        static double density = 1;
        static double miNormal = 5; //KPa
        static double kappaNormal = 6.667; //Kpa

        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa

        static double lambda0 = 1;
        int nGaussPoints = 1;
        static double initial_dp_dx = 0.0;
        static double initial_dp_dy = 0.0;
        static double initial_dp_dz = 0.0;
        static double velocityDivInitialVal = 0;
        static Dictionary<double, double[]> Solution = new Dictionary<double, double[]>();
        private static List<(INode node, IDofType dof)> watchDofs = new List<(INode node, IDofType dof)>();

        #endregion

        #region Structural model BCs and Loads

        private static List<(BC, StructuralDof[], double[][], double[])> structuralDirichletBC =
            new List<(BC, StructuralDof[], double[][], double[])>()
                {
                    (BC.BottomDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationZ }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.LeftDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.RightDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0.1,0,0}}, new double[] { 0.0 }),
                    (BC.FrontDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BC.BackDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0.1,0}}, new double[] { 0.0 }),

                };

        static double[] coordsLoad1 = new double[3] { 0.0, 0.0, 0.1 };
        static double[] coordsLoad2 = new double[3] { 0.1, 0.0, 0.1 };
        static double[] coordsLoad3 = new double[3] { 0.0, 0.1, 0.1 };
        static double[] coordsLoad4 = new double[3] { 0.1, 0.1, 0.1 };
        static double[][] loadCoords = new double[4][] { coordsLoad1, coordsLoad2, coordsLoad3, coordsLoad4 };
        static List<(BC, StructuralDof[], double[][], double[])> structuralNeumannBC =
            new List<(BC, StructuralDof[], double[][], double[])>()
            {
                (BC.TopPointFlux, new StructuralDof[1]{StructuralDof.TranslationZ}, loadCoords, new double []{1E-4 / 4d})
            };

        #endregion

        #region ToDo Orestis log task 1 
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodeDisplacementLogs = new List<(double[], string, StructuralDof, int, double[])>()
            {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        static double[] structuralMonitorNodeCoords = new double[] { 0.0, 0.0, 0.1 };
        private static int structuralMonitorID;
        static StructuralDof structuralMonitorDOF = StructuralDof.TranslationZ;

        #endregion

        #region ToDo Orestis log task 2
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpVelocityGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodeGradients",-1, new double[3][])};

        static double[] monitoredGPcoordsVelocity = new double[] { 0.004086132769345323, 0.006191771651571988, 0.09369050147682026 };

        #endregion

        #region Darcy model

        //Darcy model properties

        // original case
        //static double Sv = 7e+3; // 1/(m)
        //static double Lp = 2.7e-12; // m/(KPa sec)
        //static double LplSvl_tumor = 0; // 1/(KPa sec)
        //static double LplSvl_host = 3.75e-1; // 1/(KPa sec)
        //static double pv = 4; // kPa
        //static double pl = 0d; // KPa
        //static double k_th = 7.52e-10; // m2/(KPa sec)
        //static double k_th_tumor = 7.52e-11; // m2/(KPa sec)
        //static double k_th_host = 7.52e-13; // m2/(KPa sec)

        // simplified case
        static double Sv = 0; // 1/(m)
        static double Lp = 0; // m/(KPa sec)
        static double LplSvl_tumor = 0; // 1/(KPa sec)
        static double LplSvl_host = 0; // 1/(KPa sec)
        static double pv = 0; // kPa
        static double pl = 0d; // KPa
        static double k_th_tumor = 7.52e-6; // m2/(KPa sec)
        static double k_th_host = 7.52e-6; // m2/(KPa sec)


        #endregion

        #region Darcy BCs

        private static ConvectionDiffusionDof[] constrainedDofType = new ConvectionDiffusionDof[1] { ConvectionDiffusionDof.UnknownVariable };
        private static double[] boundaryValue = new double[1] { 0d };
        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])> pressureDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BC.BottomDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.TopDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0.1}}, boundaryValue),
                (BC.LeftDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.RightDirichlet, constrainedDofType, new double[1][]{new double[3] {0.1,0,0}}, boundaryValue),
                (BC.FrontDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BC.BackDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0.1,0}}, boundaryValue),

            };

        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])>
            pressureNeumannBC = new List<(BC, ConvectionDiffusionDof[], double[][], double[])>();

        private static double boundaryValueAllBoundaries = 0;
        #endregion

        #region Darcy Initial condition values
        // Data 1:RegionType, 2:Bcstype, 3: Region CaracteristicCoords Id, 4: Bc value
        static List<(int, int, double[][], double[])> eq78InitialConditionsList = new List<(int, int, double[][], double[])>() { (0, 0, new double[3][], new double[3]) };

        private static double initialCondition = 0; // TODO Orestis delete when obsolete

        #endregion

        #region Darcy logs

        #region ToDo Orestis log task 3 
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodePressureLogs = new List<(double[], string, StructuralDof, int, double[])>()
        {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        static double[] pressureMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.07366 };
        private static int pressureMonitorID;
        static ConvectionDiffusionDof pressureMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        #endregion

        #region ToDo Orestis log task 4
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpPressureGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodePressureGradients",-1, new double[3][])};

        static double[] monitoredGPcoordsPresGradient = new double[] { 0.004086132769345323, 0.006191771651571988, 0.09369050147682026 };

        #endregion
        #endregion

        #region Cancer Cell Density (TCell) model

        private const double dummySolidVelovity = -5;

        /// <summary>
        /// Growth rate parameter 1[mol/(m3)]
        /// </summary>
        private const double K1 = 1.74E-6; // [1/s]

        /// <summary>
        /// Growth rate parameter 2[mol/(m3)]
        /// </summary>
        private const double K2 = 8.3E-3; // [mol/(m3)]

        /// <summary>
        /// Oxygen concentration (Dependent Variable) [mol/m3]
        /// </summary>
        private const double Cox = 0.2; // [mol / m3]

        #endregion

        #region Cancer Cell Density (TCell) Boundary Conditions


        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])> tCellDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
        {(BC.TopRightBackDiriclet, constrainedDofType, new double[2][]{new double[3] {0,0,0},new double[3] {0.1,0.1,0.1}}, new double[]{500d}),};

        private static List<(BC, ConvectionDiffusionDof[], double[][], double[])>
                tCellNeumannBC = new List<(BC, ConvectionDiffusionDof[], double[][], double[])>();

        #endregion

        #region Cancer Cell Density (TCell)  Initial condition

        private double initialTCellDensity = 0d;

        #endregion

        #region Cancer Cell Density (TCell)  logs

        static List<(double[], string, StructuralDof, int, double[])> nodeTCellLogs = new List<(double[], string, StructuralDof, int, double[])>()
            {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        //static double[] tCellMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.07366 };
        static double[] tCellMonitorNodeCoords = { 0.06, 0.06, 0.06 };

        private static int tCellMonitorID;

        static ConvectionDiffusionDof tCellMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        #endregion

        public UPTdummySolution()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]

        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]

        public void MonophasicEquationModel(string fileName)
        {
            ContinuumElement3DGrowth.dT = timeStep;

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

            Dictionary<int, double[][]> velocityAtGaussPoints =
                new Dictionary<int, double[][]>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                var velocityDiv = new double[nGaussPoints][];
                for (int i1 = 0; i1 < nGaussPoints; i1++)
                {
                    velocityDiv[i1] = new double[] { 0, 0, 0 };
                };
                velocityAtGaussPoints.Add(elem.Key, velocityDiv);
            }

            Dictionary<int, double> dummyFieldCOx =
                new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                dummyFieldCOx.Add(elem.Key, Cox);
            }


            miNormal = miTumor; // TODO : remove this from here
            kappaNormal = kappaTumor;

            #region loggin (defined before model builder creation to give them nodes)

            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNodeCoords, 1e-2);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, pressureMonitorNodeCoords, 1e-2);
            tCellMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, tCellMonitorNodeCoords, 1e-2);

            var p_i = new double[(int)(totalTime / timeStep)];

            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);

            double[] gp_dut_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dvt_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dwt_dz_OverTime = new double[(int)(totalTime / timeStep)];
            var divVelocity = new List<double[]>();
            divVelocity.Add(gp_dut_dx_OverTime);
            divVelocity.Add(gp_dvt_dy_OverTime);
            divVelocity.Add(gp_dwt_dz_OverTime);


            double[] gp_div_v_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dz_Overtime = new double[(int)(totalTime / timeStep)];
            var dp_dxi = new List<double[]>();
            dp_dxi.Add(gp_dP_dx_OverTime);
            dp_dxi.Add(gp_dP_dy_OverTime);
            dp_dxi.Add(gp_dP_dz_Overtime);


            double[] tCell = new double[(int)(totalTime / timeStep)];

            int monitoredGPVelocity_elemID = -1; // TODO Orestis this will be deleeted if new logs are implemented in a right way.
            int monitoredGPpressureGrad_elemID = -1; // TODO Orestis this will be deleeted if new logs are implemented in a right way.

            //TODO Orestis: implement here one for loop for each Gp Type requested Output
            #endregion

            //Create Model For Pressure
            var pressureModel = new Eq78ModelProviderForStaggeredSolutionex7ref(comsolReader, k_th_tumor, k_th_host, Lp, Sv, pv,
                LplSvl_tumor, LplSvl_host, pl, velocityDivergenceAtElementGaussPoints, pressureMonitorID, pressureMonitorDOF, pressureDirichletBC, pressureNeumannBC);

            //Create Model For Structural
            var structuralModel = new Eq9ModelProviderForStaggeredSolutionEx7Ref(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, density, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, structuralMonitorDOF, structuralNeumannBC, structuralDirichletBC);


            //Create Model For TCell
            var tCellModel = new TCellModelProvider(K1, K2, dummyFieldCOx, velocityAtGaussPoints, comsolReader, tCellMonitorDOF, tCellMonitorID, tCellDirichletBC, tCellNeumannBC, initialTCellDensity);

            

            var equationModel = new CoupledBiphasicTCellModelProvider(pressureModel, structuralModel, tCellModel, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, velocityAtGaussPoints, timeStep,
                totalTime, incrementsPertimeStep);

            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.000000001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging

                //TODO Orestis: implement here one for loop for each "node Type" requested  log using the following commands
                monitoredGPVelocity_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoordsVelocity, 1e-1); //Todo Orestis delete these commands1
                monitoredGPpressureGrad_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoordsPresGradient, 1e-1);
                //nodal logs
                p_i[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), pressureMonitorDOF];
                //p_i[currentTimeStep] = 0d;
                //structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                //structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                //structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];
                structuralResultsX[currentTimeStep] = 0d;
                structuralResultsY[currentTimeStep] = 0d;
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), structuralMonitorDOF];

                //gp (element) logs
                gp_dP_dx_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).xcoeff_OverTimeAtGp1[0];
                gp_dP_dy_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).ycoeff_OverTimeAtGp1[0];
                gp_dP_dz_Overtime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).zcoeff_OverTimeAtGp1[0];
                gp_dut_dx_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term1[0];
                gp_dvt_dy_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term2[0];
                gp_dwt_dz_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term3[0];
                gp_div_v_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence[0];


                tCell[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[2].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[2].GetNode(tCellMonitorID), tCellMonitorDOF];

                //model maximus (DO NOT ERASE)
                //modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0])).ToArray().Max();
                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                /*if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }*/

                #endregion
                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    (equationModel.ParentAnalyzers[j] as NewmarkDynamicAnalyzer).AdvanceStep();
                }

                

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                equationModel.SaveStateFromElements();

                //Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            Assert.True(ResultChecker.CheckResults(structuralResultsZ, expectedDisplacments(), 1e-1));
            Assert.True(ResultChecker.CheckResults(p_i, expectedPressurevalues(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dut_dx_OverTime, expected_dutdx_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dvt_dy_OverTime, expected_dvtdy_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dwt_dz_OverTime, expected_dwtdz_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_div_v_OverTime, expected_div_vs_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dP_dx_OverTime, expected_dpdx_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dP_dy_OverTime, expected_dpdy_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dP_dz_Overtime, expected_dpdz_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(tCell, expected_Tc_values(), 1e-1));





            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;

            //pressure_{ pr}_F_{Fval}_e{Fval_e}_LOGGEDval_ PAth name do not erase this exampleNo_{exNo}_caseNo_{caseNo}_
            /*var writer = new MGroup.LinearAlgebra.Output.Array1DWriter();
            var patSelection = 0;
            var outputPath = patSelection == 0 ? "../../../StaggeredSolutionPresDynamex7ref/results1/" : $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\";
            writer.WriteToFile(structuralResultsX, outputPath + $@"u_{1}_.txt");
            writer.WriteToFile(structuralResultsY, outputPath + $@"u_{2}_.txt");
            writer.WriteToFile(structuralResultsZ, outputPath + $@"u_{3}_.txt");
            writer.WriteToFile(gp_dut_dx_OverTime, outputPath + $@"gp_dut_dx_OverTime.txt");
            writer.WriteToFile(gp_dvt_dy_OverTime, outputPath + $@"gp_dvt_dy_OverTime.txt");
            writer.WriteToFile(gp_dwt_dz_OverTime, outputPath + $@"gp_dwt_dz_OverTime.txt");
            writer.WriteToFile(gp_div_v_OverTime, outputPath +  $@"gp_div_v_OverTime_.txt");



            writer.WriteToFile(p_i, outputPath + $@"p_i.txt");
            writer.WriteToFile(gp_dP_dx_OverTime, outputPath + $@"gp_dP_dx_OverTime.txt");
            writer.WriteToFile(gp_dP_dy_OverTime, outputPath + $@"gp_dP_dy_OverTime.txt");
            writer.WriteToFile(gp_dP_dz_Overtime, outputPath + $@"gp_dP_dz_Overtime.txt");*/


            //var path = outputPath+"dp_dxi_mslv.csv";
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredBiphasicTCell/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredBiphasicTCell/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredBiphasicTCell/pi_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(tCell, "../../../StaggeredBiphasicTCell/tCell_nodes_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(divVelocity), "../../../StaggeredBiphasicTCell/dut_dxi_GP_mslv.csv");

        }


        public static double[] expectedDisplacments()
        {
            return new double[] {
             3.6617740000e-09,
             1.8250890000e-08,
             4.7236380000e-08,
             9.0350780000e-08,
             1.4735030000e-07,
             2.1802150000e-07,
             3.0216870000e-07,
             3.9960820000e-07,
             5.1015800000e-07,
             6.3363370000e-07,};
        }

        public static double[] expectedPressurevalues()
        {
            return new double[] {
            6.7032500000e-05,
            1.3046210000e-04,
            2.2732080000e-04,
            1.1742530000e-04,
            1.0975700000e-05,
            -9.5261250000e-05,
            -2.8098200000e-04,
            -5.0883160000e-04,
            -7.5607230000e-04,
            -1.0224200000e-03,
            };
        }

        public static double[] expected_dutdx_values()
        {
            return new double[] {
            1.5046278105920462E-05,
            8.8617497660219155E-05,
            0.00026689510562568779,
            0.00053376957178620412,
            0.00081159440063547543,
            0.0010392324512686118,
            0.0011841596351018918,
            0.0012262079260648756,
            0.0011580851478102409,
            0.00097666754328661932,
            };
        }

        public static double[] expected_dvtdy_values()
        {
            return new double[] {
            -4.9183616146580148E-05,
            -0.00031990793942135367,
            -0.0010452126289494284,
            -0.0023121523877487429,
            -0.0040322324519872065,
            -0.0061097024018783273,
            -0.0084832146516789931,
            -0.011108490603471471,
            -0.013965094795610738,
            -0.017041826768603045
            };
        }

        public static double[] expected_dwtdz_values()
        {
            return new double[] {
            0.024013273060017812,
            0.095527520745864361,
            0.18931185141078608,
            0.28069438945206548,
            0.36992113985726005,
            0.4572951803247241,
            0.54300618640408982,
            0.62718236463256549,
            0.70985946584687643,
            0.79103569439181676
            };
        }

        public static double[] expected_div_vs_values()
        {
            return new double[] {
            0.023979135721977154,
            0.095296230304103224,
            0.18853353388746233,
            0.27891600663610294,
            0.3667005018059083,
            0.45222471037411438,
            0.53570713138751269,
            0.61730008195515884,
            0.697052456199076,
            0.77497053516650038
            };
        }

        public static double[] expected_dpdx_values()
        {
            return new double[] {
           -0.51071856532594151,
            -1.7887001011389974,
            -3.4085533329685966,
            -4.3925950633681392,
            -5.2590821460816093,
            -5.8433941491425614,
            -6.2892322572273516,
            -6.6400814219331838,
            -6.91972551191365,
            -7.1456532816455649
            };
        }

        public static double[] expected_dpdy_values()
        {
            return new double[10];
        }

        public static double[] expected_dpdz_values()
        {
            return new double[10];
        }

        public static double[] expected_Tc_values()
        {
            return new double[] {
            0.39161779623154619,
            0.77197862637823178,
            1.14114595438075,
            1.4991842343623834,
            1.8461589027996081,
            2.1821363705976768,
            2.507184015072121,
            2.8213701718371418,
            3.1247641266018533,
            3.4174361068753587,
            };
        }







    }
}
