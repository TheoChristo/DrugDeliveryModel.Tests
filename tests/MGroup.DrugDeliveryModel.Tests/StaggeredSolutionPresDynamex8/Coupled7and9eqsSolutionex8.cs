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

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled7and9eqsSolutionex8
    {
        const double Sc = 0.1;

        private const double timeStep = 0.00001; // in sec
        const double totalTime = 0.00040;//0.00100; // in sec
        static int incrementsPertimeStep = 1;


        // strucutral model Loads
        static StructuralDof loadedDof = StructuralDof.TranslationZ;
        //static double load_value = 0.0001 / 4; //[kN]
        static double load_value = 0.0001 / 4; //[kN]


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


        private static int pressureMonitorID;

        static double divPMonitorGPX = 0.01714618146;
        static double divPMonitorGPY = 0.05114110504;
        static double divPMonitorGPZ = 0.06533804930;

        //static double divPMonitorGPX = 0.05;
        //static double divPMonitorGPY = 0.05;
        //static double divPMonitorGPZ = 0.05;

        //static double divPMonitorGPX = 0;
        //static double divPMonitorGPY = 0;
        //static double divPMonitorGPZ = 0.1;

        //static double divPMonitorGPX = 0.004086132769345323;
        //static double divPMonitorGPY = 0.006191771651571988;
        //static double divPMonitorGPZ = 0.09369050147682026;

        //GP with coordinates: 0.05301208792514899 0.053825572057669926 0.052065045951539365 is in element with id: 1251 NEAR [0.5, 0.5, 0.5]
        //GP with coordinates: 0.004086132769345323 0.006191771651571988 0.09369050147682026 is in element with id: 740  NEAR [0,0,0.1]


        //static double[] monitoredGPcoords = new double[] { 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 };
        static double[] monitoredGPcoords = new double[] { 0.004086132769345323, 0.006191771651571988, 0.09369050147682026 };

        private static int divPMonitorID;

        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationZ;

        //Darcy model properties

        // original case
        //static double Sv = 7e+3; // 1/(m)
        //static double Lp = 2.7e-12; // m/(KPa sec)
        //static double LplSvl = 0; // 1/(KPa sec)
        //static double pv = 4; // kPa
        //static double pl = 0d; // KPa
        //static double k_th = 7.52e-10; // m2/(KPa sec)

        // simplified case
        static double Sv = 0; // 1/(m)
        static double Lp = 0; // m/(KPa sec)
        static double LplSvl = 0; // 1/(KPa sec)
        static double pv = 0; // kPa
        static double pl = 0d; // KPa
        static double k_th = 7.52e-2; // m2/(KPa sec)




        //static int nodeIdToMonitor = 36;

        //The above coordinates are in the reference configuration.
        //The actual coordinates will be computed in the analysis.

        //static double structuralMonitorNodeX = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        //static double structuralMonitorNodeY = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa
        //static double structuralMonitorNodeZ = 0.05;//0.05; einai gia artio arithmo diakritopoishshs me hexa

        static double structuralMonitorNodeX = 0.0;
        static double structuralMonitorNodeY = 0.0;
        static double structuralMonitorNodeZ = 0.1;

        //static double[] structuralMonitorNodeCoords = new double[]
        //    { 0.04930793848882013,0.04994681648346263,0.04953188199244812 };
        //static double[] structuralMonitorNodeCoords = new double[]
        //    { 0.04379781548934232, 0.0492383785341323, 0.02519415945383001 };
        static double[] structuralMonitorNodeCoords = new double[]
            { 0.0, 0.0, 0.1 };

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

        static double[] expectedP_GP = new double[] { -6.256436420041747e-07, -1.3547937739021014e-06, -2.134649449582765e-06, -2.938416080376368e-06, -3.752297639083258e-06, -4.568996435665539e-06, -5.384465083785849e-06, -6.196283763502672e-06, -7.002849729927342e-06, -7.802972740849993e-06, -8.595673338476355e-06, -9.380082501804736e-06, -1.0155391952036858e-05, -1.0920829763956986e-05, -1.1675648615292048e-05, -1.2419120342525798e-05, -1.3150533638329165e-05, -1.3869193308369004e-05, -1.457442029617367e-05, -1.526555208004013e-05, -1.594194324351657e-05, -1.660296611976671e-05, -1.7248011459493847e-05, -1.787648909680102e-05, -1.8487828599733248e-05, -1.908147989845301e-05, -1.965691388712747e-05, -2.021362299719242e-05, -2.075112174048476e-05, -2.126894722116352e-05, -2.176665961559485e-05, -2.2243842619525905e-05, -2.2700103861972035e-05, -2.3135075285330776e-05, -2.354841349130354e-05, -2.393980005226653e-05, -2.4308941787803043e-05, -2.4655571006159904e-05, -2.497944571046521e-05, -2.5280349769591517e-05, -2.555809305361652e-05, -2.581251153389252e-05, -2.6043467347789186e-05, -2.6250848828245337e-05, -2.643457049831271e-05, -2.65945730309362e-05, -2.673082317427547e-05, -2.68433136429236e-05, -2.693206297543751e-05, -2.6997115358644027e-05, -2.7038540419248e-05, -2.705643298330093e-05, -2.705091280416254e-05, -2.7022124259615675e-05, -2.6970236018849343e-05, -2.6895440680066232e-05, -2.6797954379520624e-05, -2.667801637282005e-05, -2.6535888589372265e-05, -2.637185516089257e-05, -2.618622192491657e-05, -2.597931590429914e-05, -2.5751484763708963e-05, -2.5503096244149318e-05, -2.523453757656614e-05, -2.4946214875623772e-05, -2.4638552514744853e-05, -2.431199248353232e-05, -2.396699372870066e-05, -2.360403147965986e-05, -2.322359655989891e-05, -2.282619468532536e-05, -2.2412345750724504e-05, -2.198258310549241e-05, -2.1537452819812016e-05, -2.107751294242072e-05, -2.0603332751128652e-05, -2.011549199722591e-05, -1.9614580144916535e-05, -1.910119560689576e-05, -1.8575944977185543e-05, -1.8039442262312528e-05, -1.7492308111907e-05, -1.693516904977495e-05, -1.6368656706476152e-05, -1.5793407054421053e-05, -1.5210059646464891e-05, -1.4619256858970884e-05, -1.4021643140264752e-05, -1.3417864265392125e-05, -1.280856659805457e-05, -1.2194396360571182e-05, -1.1575998912680843e-05, -1.095401803996935e-05, -1.0329095252678816e-05, -9.701869095614402e-06, -9.07297446983909e-06, -8.443041966814352e-06, -7.812697215601512e-06 };
        static double[] expectedDisplacementZNode = new double[] { 7.363668893030732e-09, 2.5761887877319318e-08, 5.70115014874203e-08, 1.0199417962869543e-07, 1.6111562324654824e-07, 2.345350934186782e-07, 3.222801854819649e-07, 4.243042747939722e-07, 5.405153302312236e-07, 6.707904347588978e-07, 8.149831782450316e-07, 9.72927502706326e-07, 1.144439788712001e-06, 1.3293200764438096e-06, 1.5273528675335897e-06, 1.7383077302170753e-06, 1.9619398185917123e-06, 2.1979903608989367e-06, 2.446187143810359e-06, 2.706245005712201e-06, 2.97786634497212e-06, 3.260741645648848e-06, 3.554550021324733e-06, 3.858959776833323e-06, 4.173628987183327e-06, 4.4982060927292285e-06, 4.832330509498858e-06, 5.175633253505171e-06, 5.527737577816046e-06, 5.888259621119012e-06, 6.256809066490117e-06, 6.632989809054993e-06, 7.016400631213204e-06, 7.40663588408404e-06, 7.803286173822064e-06, 8.205939051444086e-06, 8.614179704805687e-06, 9.027591651364676e-06, 9.44575743037129e-06, 9.86825929313e-06, 1.0294679889986136e-05, 1.0724602952701256e-05, 1.1157613970895282e-05, 1.1593300861249765e-05, 1.203125462818632e-05, 1.2471070014756037e-05, 1.291234614250065e-05, 1.3354687139073446e-05, 1.3797702752437716e-05, 1.4241008950492924e-05, 1.4684228505013333e-05, 1.512699155882082e-05, 1.5568936175152747e-05, 1.6009708868226998e-05, 1.6448965114049556e-05, 1.6886369840555107e-05, 1.732159789621823e-05, 1.7754334496321283e-05, 1.8184275646115336e-05, 1.861112854016225e-05, 1.9034611937199034e-05, 1.945445650991976e-05, 1.9870405169125824e-05, 2.0282213361751447e-05, 2.0689649342328735e-05, 2.1092494417514197e-05, 2.1490543163357156e-05, 2.1883603615049627e-05, 2.227149742895581e-05, 2.265406001677948e-05, 2.3031140651786245e-05, 2.3402602547057626e-05, 2.3768322905812383e-05, 2.4128192943889832e-05, 2.448211788454774e-05, 2.483001692578529e-05, 2.5171823180458423e-05, 2.5507483589510936e-05, 2.5836958808699867e-05, 2.616022306924752e-05, 2.647726401290536e-05, 2.6788082501966378e-05, 2.7092692404812564e-05, 2.739112035763236e-05, 2.7683405502990065e-05, 2.796959920597384e-05, 2.824976474869246e-05, 2.8523977003931938e-05, 2.8792322088822693e-05, 2.9054896999404773e-05, 2.9311809227013853e-05, 2.9563176357443383e-05, 2.9809125653868805e-05, 3.004979362454759e-05, 3.028532557633501e-05, 3.051587515507827e-05, 3.0741603873972754e-05, 3.09626806309822e-05, 3.117928121644024e-05 };

        public Coupled7and9eqsSolutionex8()
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
        //[InlineData("../../../DataFiles/workingQHexaMesh6x6x6_1Domain.mphtxt")]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
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

            #region loggin (defined before model builder creation to give them nodes)
            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNodeCoords, 1e-2);
            //pressureMonitorID = Utilities.FindRandomInternalNode(comsolReader.NodesDictionary, modelMinX, modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new double[] { 0.055, 0.0559, 0.07366 }, 1e-2);
            //pressureMonitorID = structuralMonitorID;

            var p_i = new double[(int)(totalTime / timeStep)];

            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);

            //double[] gp_dut_dx_OverTime = new double[(int)(totalTime / timeStep)];
            //double[] gp_dvt_dy_OverTime = new double[(int)(totalTime / timeStep)];
            //double[] gp_dwt_dz_OverTime = new double[(int)(totalTime / timeStep)];
            //var divVelocity = new List<double[]>();
            //divVelocity.Add(gp_dut_dx_OverTime);
            //divVelocity.Add(gp_dvt_dy_OverTime);
            //divVelocity.Add(gp_dwt_dz_OverTime);


            double[] gp_div_v_OverTime = new double[(int)(totalTime / timeStep)];
            //double[] gp_dP_dx_OverTime = new double[(int)(totalTime / timeStep)];
            //double[] gp_dP_dy_OverTime = new double[(int)(totalTime / timeStep)];
            //double[] gp_dP_dz_Overtime = new double[(int)(totalTime / timeStep)];
            //var dp_dxi = new List<double[]>();
            //dp_dxi.Add(gp_dP_dx_OverTime);
            //dp_dxi.Add(gp_dP_dy_OverTime);
            //dp_dxi.Add(gp_dP_dz_Overtime);

            int monitoredGP_elemID = -1;
            #endregion




            var eq78Model = new Eq78ModelProviderForStaggeredSolutionex8(comsolReader, k_th, Lp, Sv, pv, LplSvl, pl,
                velocityDivergenceAtElementGaussPoints,

                boundaryValueAllBoundaries, initialCondition, pressureMonitorID, eq7n8dofTypeToMonitor, modelMinX,
                modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            //COMMITED BY NACHO
            //jkkk bn///////vji typ[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[00u-----------------------------------
            var eq9Model = new Eq9ModelProviderForStaggeredSolutionex8(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, eq9dofTypeToMonitor, loadedDof, load_value, modelMinX, modelMaxX, modelMinY,
                modelMaxY, modelMinZ, modelMaxZ);

            var equationModel = new Coupled7and9eqsModelex8(eq78Model, eq9Model, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, timeStep,
                totalTime, incrementsPertimeStep);






            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging
                monitoredGP_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoords, 1e-1);

                //nodal logs
                p_i[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), ConvectionDiffusionDof.UnknownVariable];
                //p_i[currentTimeStep] = 0d;
                //structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                //structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                //structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];
                structuralResultsX[currentTimeStep] = 0d;
                structuralResultsY[currentTimeStep] = 0d;
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];

                //gp (element) logs
                //gp_dP_dx_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).xcoeff_OverTimeAtGp1[0];
                //gp_dP_dy_OverTime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).ycoeff_OverTimeAtGp1[0];
                //gp_dP_dz_Overtime[currentTimeStep] = ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGP_elemID]).zcoeff_OverTimeAtGp1[0];
                //gp_dut_dx_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term1[0];
                //gp_dvt_dy_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term2[0];
                //gp_dwt_dz_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence_term3[0];
                //gp_div_v_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGP_elemID]).velocityDivergence[0];

                //model maximus (DO NOT ERASE)
                //modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0])).ToArray().Max();
                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints.Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                if (Solution.ContainsKey(currentTimeStep))
                {
                    // Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    // Solution.Add(currentTimeStep, allValues);
                }

                #endregion

                (equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                //equationModel.SaveStateFromElements();

                //Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;

            //pressure_{ pr}_F_{Fval}_e{Fval_e}_LOGGEDval_ PAth name do not erase this exampleNo_{exNo}_caseNo_{caseNo}_
            var writer = new MGroup.LinearAlgebra.Output.Array1DWriter();
            var patSelection = 0;
            var outputPath = patSelection == 0 ? "../../../StaggeredSolutionPresDynamex8/results1/" : $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\";
            writer.WriteToFile(structuralResultsX, outputPath + $@"u_{1}_.txt");
            writer.WriteToFile(structuralResultsY, outputPath + $@"u_{2}_.txt");
            writer.WriteToFile(structuralResultsZ, outputPath + $@"u_{3}_.txt");
            //writer.WriteToFile(gp_dut_dx_OverTime, outputPath + $@"gp_dut_dx_OverTime.txt");
            //writer.WriteToFile(gp_dvt_dy_OverTime, outputPath + $@"gp_dvt_dy_OverTime.txt");
            //writer.WriteToFile(gp_dwt_dz_OverTime, outputPath + $@"gp_dwt_dz_OverTime.txt");
            //writer.WriteToFile(gp_div_v_OverTime, outputPath +  $@"gp_div_v_OverTime_.txt");



            writer.WriteToFile(p_i, outputPath + $@"p_i.txt");
            //writer.WriteToFile(gp_dP_dx_OverTime, outputPath + $@"gp_dP_dx_OverTime.txt");
            //writer.WriteToFile(gp_dP_dy_OverTime, outputPath + $@"gp_dP_dy_OverTime.txt");
            //writer.WriteToFile(gp_dP_dz_Overtime, outputPath + $@"gp_dP_dz_Overtime.txt");


            //var path = outputPath+"dp_dxi_mslv.csv";
            //CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredSolutionPresDynamex8/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredSolutionPresDynamex8/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredSolutionPresDynamex8/pi_nodes_mslv.csv");
            //CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(divVelocity), "../../../StaggeredSolutionPresDynamex8/dut_dxi_GP_mslv.csv");

            Assert.True(CompareResults(p_i, displacements));
        }

        private bool CompareResults(double[] pressure, List<double[]> displacements)
        {
            bool ret = true;
            double tolerance = 1e-4;

            for (var i = 0; i < pressure.Length; i++)
            {
                double pError = Math.Abs(pressure[i] - expectedP_GP[i]);
                double dispZError = Math.Abs(displacements[2][i] - expectedDisplacementZNode[i]);

                if (pError > tolerance || dispZError > tolerance)
                {
                    Console.WriteLine("Wrong result on step " + i);
                    Console.WriteLine("pError: " + pError);
                    Console.WriteLine("dispZError " + dispZError);
                    ret = false;
                    break;
                }
            }

            return ret;
        }




    }
}
