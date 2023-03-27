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
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;
using TriangleNet;
using TriangleNet.Meshing;
using Castle.DynamicProxy.Generators.Emitters.SimpleAST;
using System.Reflection.PortableExecutable;
using static Xunit.Assert;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled78_9_13VanillaSourceSolution
    {
        const double Sc = 0.1;

        private const double timeStep = 1E-5; // in sec
        const double totalTime = 3E-4; // in sec
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

        //Array of doublle arrys that contain various node coordintates (in the same face)
        private double[][] test = new double[2][] { new double[] { 1, 1, 1 }, new double[] { 1, 1, 1 } };
        
        #region Structural model BCs and Loads
        
        private static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralDirichletBC = 
            new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])>()
                {
                    (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BottomDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationZ }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.LeftDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationX }, new double[1][]{new double[3] {0.1,0,0}}, new double[] { 0.0 }),
                    (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.FrontDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0,0}}, new double[] { 0.0 }),
                    (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BackDirichlet,
                        new StructuralDof[1] { StructuralDof.TranslationY }, new double[1][]{new double[3] {0,0.1,0}}, new double[] { 0.0 }),
                    
                };
        
        static double[] coordsLoad1 = new double[3] { 0.0, 0.0, 0.1 };
        static double[] coordsLoad2 = new double[3] { 0.1, 0.0, 0.1 };
        static double[] coordsLoad3 = new double[3] { 0.0, 0.1, 0.1 };
        static double[] coordsLoad4 = new double[3] { 0.1, 0.1, 0.1 };
        static double[][] loadCoords = new double[4][] { coordsLoad1, coordsLoad2, coordsLoad3, coordsLoad4 };
        static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])> structuralNeumannBC =
            new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, StructuralDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopPointFlux, new StructuralDof[1]{StructuralDof.TranslationZ}, loadCoords, new double []{-1E-4 / 4d})
            };
        
        #endregion

        #region Solid Displacement Logging
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodeDisplacementLogs = new List<(double[], string, StructuralDof, int, double[])>()
        {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};

        static double[] structuralMonitorNodeCoords = new double[]
            { 0.0, 0.0, 0.1 };
        private static int structuralMonitorID;
        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationZ;
        #endregion

        #region Solid Velocity Logging
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpVelocityGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodeGradients",-1, new double[3][])};

        //10 step  u - P coupling
        //static double[] monitoredGPCoordsVelocity = new double[] { 0.09, 0.09, 0.09 };
        
        //30 step - fluid velocity test
        static double[] monitoredGPCoordsVelocity = new double[] { 0.08, 0.08, 0.08 };

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
        private static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> pressureDirichletBC = 
            new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BottomDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0.1}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.LeftDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet, constrainedDofType, new double[1][]{new double[3] {0.1,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.FrontDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0,0}}, boundaryValue),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BackDirichlet, constrainedDofType, new double[1][]{new double[3] {0,0.1,0}}, boundaryValue),
                    
            };

        private static List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>
            pressureNeumannBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double [])>();
        
        private static double boundaryValueAllBoundaries = 0;
        #endregion

        #region Darcy Initial condition values
        // Data 1:RegionType, 2:Bcstype, 3: Region CaracteristicCoords Id, 4: Bc value
        static List<(int, int, double[][], double[])> eq78InitialConditionsList = new List<(int, int, double[][], double[])>() {(0, 0, new double[3][], new double[3])};

        private static double initialCondition = 0; // TODO Orestis delete when obsolete
        
        #endregion

        #region Darcy logs

        #region Fluid Pressure Logging 
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: logged Doftype, 4: found node Id, 5: results 
        /// </summary>
        static List<(double[], string, StructuralDof, int, double[])> nodePressureLogs = new List<(double[], string, StructuralDof, int, double[])>()
        {(new double[]{ 0.04930793848882013,0.04994681648346263,0.075 }, "CornerNodeTranslationZ.txt",StructuralDof.TranslationZ,-1, new double[0])};
        
        //10 step  u - P coupling
        static double[] pressureMonitorNodeCoords = new double[] { 0.055, 0.0559, 0.05 };
        
        private static int pressureMonitorID;
        
        static ConvectionDiffusionDof eq7n8dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;
        
        #endregion

        #region Fluid Pressure Gradient & Velocity Logging
        /// <summary>
        /// Data 1:coords, 2:outputfileString, 3: found element Id, 4: results 
        /// </summary>
        static List<(double[], string, int, double[][])> gpPressureGraadientLogs = new List<(double[], string, int, double[][])>()
        {(new double[]{ 0.05301208792514899, 0.053825572057669926, 0.052065045951539365 }, "CenterNodePressureGradients",-1, new double[3][])};

        //10 step  u - P coupling
        //static double[] monitoredGPcoordsPresGradient = new double[] { 0.055, 0.0559, 0.05 };
        
        //30 step - fluid velocity test
        static double[] monitoredGPcoordsPresGradient = new double[] { 0.08, 0.08, 0.08 };
        
        //30 step - fluid velocity test
        static double[] monitoredGPcoordsFluidVelocity = new double[] { 0.08, 0.08, 0.08 };
        static int fluidVelocityMonitorID;

        #endregion
        #endregion

        #region CoxParams
        /// <summary>
        /// The average value of the three components of the fluid velocity vector  [m/s]
        /// </summary>
        private Dictionary<int, double[]> FluidSpeed = new Dictionary<int, double[]>(); // 2.32E-4 [m/s]
        const double FluidSpeedInit = 0;//2.32E-4;

        static double SvCox = 7E3; // 1 / m
        /// <summary>
        /// Diffusivity of oxygen [m2/s]
        /// </summary>
        private const double Dox = 1.79E-5; // [m2/s]

        /// <summary>
        /// Oxygen uptake [mol/(m3*s)]
        /// </summary>
        private const double Aox = 2.55E-2; // [mol/(m3*s)]

        /// <summary>
        /// Oxygen uptake [mol/m3]
        /// </summary>
        private const double Kox = 4.64E-3; // [mol / m3]

        /// <summary>
        /// Oxygen permeability across tumor vessel walls [m/s]
        /// </summary>
        private const double PerOx = 3.55E-4; // [m/s]

        /// <summary>
        /// Initial Oxygen Concentration [mol/m3]
        /// </summary>
        private const double CInitialOx = 0d; // [mol/m3]

        private const double CInitOx =0.2; // [mol/m3]
        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private Dictionary<int, double> T = new Dictionary<int, double>();// 500 [cells]
        const double TInit = 500;

        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof coxMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        private double[] coxMonitorNodeCoords = { 0.09, 0.0, 0.05 };
        private double[] vfMonitorGpCoords = { 0.09, 0.09, 0.09 };
        private static int coxMonitorID;
        private static int vfMonitorGpID;
        private List<double[]> vf_calculated = new List<double[]>();

        private readonly Func<double> independentLinearSource = () => PerOx * SvCox * CInitOx;
        
        private readonly Func<double> dependentLinearSource = () => -PerOx * SvCox;
        

        #endregion

        #region CoxBCs
        List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionDirichletBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.TopDirichlet, new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2}),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.RightDirichlet,new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2d}),
                (BoundaryAndInitialConditionsUtility.BoundaryConditionCase.BackDirichlet, new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[1][]{new double[3] {0.1, 0.1, 0.1}}, new double[] {0.2d}),
            };
        List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> convectionDiffusionNeumannBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();
        #endregion


        public Coupled78_9_13VanillaSourceSolution()
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
            ContinuumElement3DGrowth.dT = timeStep;

            //Read geometry
            var comsolReader = new ComsolMeshReader(fileName);

            // initialize Shared quantities of Coupled model
            Dictionary<int, double> lambda = new Dictionary<int, double>(comsolReader.ElementConnectivity.Count());
            foreach (var elem in comsolReader.ElementConnectivity)
            {
                lambda.Add(elem.Key, lambda0);
                
                //Cox Init
                FluidSpeed.Add(elem.Key, new double[] { FluidSpeedInit, FluidSpeedInit, FluidSpeedInit });
                T.Add(elem.Key, TInit);
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

            miNormal = miTumor; // TODO : remove this from here
            kappaNormal = kappaTumor;

            #region loggin (defined before model builder creation to give them nodes)
            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNodeCoords, 1e-2);
            //pressureMonitorID = Utilities.FindRandomInternalNode(comsolReader.NodesDictionary, modelMinX, modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, pressureMonitorNodeCoords , 1e-2);
            //pressureMonitorID = structuralMonitorID;
            coxMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, coxMonitorNodeCoords, 1e-2);
            fluidVelocityMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, monitoredGPcoordsFluidVelocity, 1e-2);
            var p_i = new double[(int)(totalTime / timeStep)];

            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>
            {
                structuralResultsX,
                structuralResultsY,
                structuralResultsZ
            };

            double[] solidVelocityResultsX = new double[(int)(totalTime / timeStep)];
            double[] solidVelocityResultsY = new double[(int)(totalTime / timeStep)];
            double[] solidVelocityResultsZ = new double[(int)(totalTime / timeStep)];
            var solidVelocities = new List<double[]>();
            solidVelocities.Add(solidVelocityResultsX);
            solidVelocities.Add(solidVelocityResultsY);
            solidVelocities.Add(solidVelocityResultsZ);

            double[] gp_dut_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dvt_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dwt_dz_OverTime = new double[(int)(totalTime / timeStep)];
            var divVelocity = new List<double[]>
            {
                gp_dut_dx_OverTime,
                gp_dvt_dy_OverTime,
                gp_dwt_dz_OverTime
            };
            
            double[] gp_div_v_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dx_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dy_OverTime = new double[(int)(totalTime / timeStep)];
            double[] gp_dP_dz_Overtime = new double[(int)(totalTime / timeStep)];
            var dp_dxi = new List<double[]>();
            dp_dxi.Add(gp_dP_dx_OverTime);
            dp_dxi.Add(gp_dP_dy_OverTime);
            dp_dxi.Add(gp_dP_dz_Overtime);
            
            double[] uFluid_t = new double[(int)(totalTime / timeStep)];
            double[] vFluid_t = new double[(int)(totalTime / timeStep)];
            double[] wFluid_t = new double[(int)(totalTime / timeStep)];
            var fluidVelocity = new List<double[]>();
            fluidVelocity.Add(uFluid_t);
            fluidVelocity.Add(vFluid_t);
            fluidVelocity.Add(wFluid_t);
            
            double[] coxResults = new double[(int)(totalTime / timeStep)];

            int monitoredGPVelocity_elemID = -1; // TODO Orestis this will be deleeted if new logs are implemented in a right way.
            int monitoredGPpressureGrad_elemID = -1; // TODO Orestis this will be deleeted if new logs are implemented in a right way.

            //TODO Orestis: implement here one for loop for each Gp Type requested Output
            #endregion


            //Create Model For Pressure
            var eq78Model = new Eq78ModelProviderForStaggeredSolutionex7ref(comsolReader, k_th_tumor, k_th_host, Lp, Sv, pv,
                LplSvl_tumor, LplSvl_host, pl, velocityDivergenceAtElementGaussPoints, pressureMonitorID, eq7n8dofTypeToMonitor,pressureDirichletBC, pressureNeumannBC);

            //Create Model For Structural
            var eq9Model = new Eq9ModelProviderForStaggeredSolutionEx7Ref(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, density, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, eq9dofTypeToMonitor, structuralNeumannBC, structuralDirichletBC);

            //Create Model For Oxygen
            var coxModel = new CoxVanillaSourceModelBuilder(comsolReader, FluidSpeed, independentLinearSource, dependentLinearSource, Dox, Aox, Kox, PerOx, SvCox, CInitialOx, 
                                            
                                            coxMonitorID, coxMonitorDOF, convectionDiffusionDirichletBC, convectionDiffusionNeumannBC);
           

            //COMMITED BY NACHO 
            //jkkk bn///////vji typ[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[00u-----------------------------------

            var equationModel = new Coupled78_9_13_VanillaSourceModel(eq78Model, coxModel, eq9Model , comsolReader, lambda, pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, FluidSpeed, k_th_tumor, timeStep,
                totalTime, incrementsPertimeStep);
            
            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers, equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance:1E-5);                                                       
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging
                
                //TODO Orestis: implement here one for loop for each "node Type" requested  log using the following commands
                monitoredGPVelocity_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoordsFluidVelocity, 1e-1); //Todo Orestis delete these commands1
                monitoredGPpressureGrad_elemID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], monitoredGPcoordsPresGradient, 1e-1);
                vfMonitorGpID = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[2], vfMonitorGpCoords, 1e-1);

                //nodal logs
                p_i[currentTimeStep] =((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), eq7n8dofTypeToMonitor];
                //p_i[currentTimeStep] = 0d;
                //structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                //structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                //structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];
                structuralResultsX[currentTimeStep] = 0d;
                structuralResultsY[currentTimeStep] = 0d;
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), eq9dofTypeToMonitor];

                //gp (element) logs
                gp_dP_dx_OverTime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).xcoeff_OverTimeAtGp1[0];
                gp_dP_dy_OverTime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).ycoeff_OverTimeAtGp1[0];
                gp_dP_dz_Overtime[currentTimeStep] =((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[monitoredGPpressureGrad_elemID]).zcoeff_OverTimeAtGp1[0];
                gp_dut_dx_OverTime[currentTimeStep]= ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term1[0];
                gp_dvt_dy_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term2[0];
                gp_dwt_dz_OverTime[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence_term3[0];
                gp_div_v_OverTime[currentTimeStep]= ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocityDivergence[0];

                var fluidVelocityX = coxModel.FluidSpeed[monitoredGPVelocity_elemID][0];
                var fluidVelocityY = coxModel.FluidSpeed[monitoredGPVelocity_elemID][1];
                var fluidVelocityZ = coxModel.FluidSpeed[monitoredGPVelocity_elemID][2];
                uFluid_t[currentTimeStep] = fluidVelocityX;
                vFluid_t[currentTimeStep] = fluidVelocityY;
                wFluid_t[currentTimeStep] = fluidVelocityZ;
                
                solidVelocityResultsX[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][0];
                solidVelocityResultsY[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][1];
                solidVelocityResultsZ[currentTimeStep] = ((ContinuumElement3DGrowth)equationModel.model[1].ElementsDictionary[monitoredGPVelocity_elemID]).velocity[0][2];
                
                
                coxResults[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[2].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[2].GetNode(coxMonitorID), coxMonitorDOF];
                //vf_calculated.Add(FluidSpeed[vfMonitorGpID]);
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

                (equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[2] as NewmarkDynamicAnalyzer).AdvanceStep();

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                equationModel.SaveStateFromElements();

                //Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../Coupling78_9_13/results/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../Coupling78_9_13/results/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../Coupling78_9_13/results/pi_nodes_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(divVelocity), "../../../Coupling78_9_13/results/dut_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(new List<double[]>() { coxResults}), "../../../Coupling78_9_13/results/cox_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(fluidVelocity), "../../../Coupling78_9_13/results/vFluid_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(solidVelocities), "../../../Coupling78_9_13/results/vSolid_GP_mslv.csv");

            
            Assert.True(ResultChecker.CheckResults(structuralResultsZ, expectedDisplacments(), 1e-3));
            
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsX, expectedSolidVelocityX(), 1e-3));
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsY, expectedSolidVelocityY(), 1e-3));
            Assert.True(ResultChecker.CheckResults(solidVelocityResultsZ, expectedSolidVelocityZ(), 1e-3));
            
            Assert.True(ResultChecker.CheckResults(p_i, expectedPressurevalues(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dut_dx_OverTime, expected_dutdx_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dvt_dy_OverTime, expected_dvtdy_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dwt_dz_OverTime, expected_dwtdz_values(), 1e-3));
            //Assert.True(ResultChecker.CheckResults(gp_div_v_OverTime, expected_div_vs_values(), 1e-1));
            Assert.True(ResultChecker.CheckResults(gp_dP_dx_OverTime, expected_dpdx_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dP_dy_OverTime, expected_dpdy_values(), 1e-3));
            Assert.True(ResultChecker.CheckResults(gp_dP_dz_Overtime, expected_dpdz_values(), 1e-3));
            //Assert.True(ResultChecker.CheckResults(coxResults, expectedCox(), 1E-1));


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

        }

        static double[] expectedCox()
        {
            return new double[]
            {
                0.0,
                1.170204359186622e-08,
                3.4570848140523114e-08,
                6.225601204495958e-08,
                8.608086597137913e-08,
                9.686224086772193e-08,
                8.581659700450394e-08,
                4.497588499747337e-08,
                -3.265279028121877e-08,
                -1.5305625845604175e-07,
/*                -3.212450237753937e-07,
                -5.413081587338103e-07,
                -8.164754991650073e-07,
                -1.1491786627453075e-06,
                -1.5411072709959789e-06,
                -1.9932591861453542e-06,
                -2.5059847270519494e-06,
                -3.0790253258415123e-06,
                -3.711547261089456e-06,
                -4.402171127147807e-06,
                -5.148997659910326e-06,
                -5.949630478759389e-06,
                -6.801196241368415e-06,
                -7.700362649817535e-06,
                -8.643354694856509e-06,
                -9.625969480152413e-06,
                -1.0643589928909174e-05,
                -1.16911976403368e-05,
                -1.2763385132388356e-05,
                -1.3854367679219798e-05,
                -1.4957994926515588e-05,
                -1.6067762444996685e-05,
                -1.7176823361648006e-05,
                -1.8278000189476375e-05,
                -1.936379695959752e-05,
                -2.04264117442482e-05,
                -2.1457749645542068e-05,
                -2.244943631247986e-05,
                -2.339283203766632e-05,
                -2.4279046475297077e-05,
                -2.5098954013120333e-05,
                -2.5843209823155e-05,
                -2.6502266608803707e-05,
                -2.7066392059629043e-05,
                -2.7525687019297424e-05,
                -2.787010436699178e-05,
                -2.8089468607896937e-05,
                -2.8173496164078403e-05,
                -2.8111816353174065e-05,
                -2.789399303873408e-05,
                -2.75095469327707e-05,
                -2.6947978528065367e-05,
                -2.6198791634939105e-05,
                -2.5251517494633797e-05,
                -2.409573943898603e-05,
                -2.272111806384629e-05,
                -2.1117416881579355e-05, 
                -1.9274528415990708e-05, 
                -1.7182500701181277e-05, 
                -1.4831564144116403e-05, 
                -1.2212158709044051e-05, 
                -9.314961380403575e-06, 
                -6.130913859454398e-06, 
                -2.6512504485408324e-06, 
                1.132473924316709e-06, 
                5.228355588917331e-06, 
                9.64411397389559e-06, 
                1.4387063455639645e-05, 
                1.946408526631728e-05, 
                2.4881599512947062e-05, 
                3.064553736013198e-05, 
                3.676131342961979e-05, 
                4.323379847030045e-05, 
                5.006729235264707e-05, 
                5.72654974418942e-05, 
                6.483149240442061e-05, 
                7.276770650191279e-05, 
                8.107589442789177e-05, 
                8.975711174111907e-05, 
                9.88116909502036e-05, 
                0.00010823921830349637, 
                0.00011803851133800888, 
                0.000128207597240658, 
                0.00013874369207462242, 
                0.00014964318092300927,
                0.0001609015990013394, 
                0.00017251361378959422, 
                0.00018447300823374915, 
                0.00019677266506577064, 
                0.00020940455229010067, 
                0.00022235970988354397, 
                0.00023562823775436208, 
                0.00024919928500516126, 
                0.0002630610405428912, 
                0.0002772007250778919, 
                0.0002916045845525636, 
                0.00030625788503870665, 
                0.0003211449091410782, 
                0.000336248953943072, 
                0.00035155233052878125*/
            };
        }

        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //Node : 0.0, 0.0, 0.1
        public static double[] expectedDisplacments()
        {
            return new double[] {
                -3.66177E-09,
                -1.82506E-08,
                -4.72353E-08,
                -9.03487E-08,
                -1.47347E-07,
                -2.18019E-07,
                -3.02169E-07,
                -3.99613E-07,
                -5.10167E-07,
                -6.33648E-07,
            }; 
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityX()
        {
            return new double[] {
                -2.57948E-08,
                -1.5538E-07 ,
                -4.79338E-07,
                -1.00336E-06,
                -1.65214E-06,
                -2.38495E-06,
                -3.20718E-06,
                -4.12942E-06,
                -5.16614E-06,
                -6.33261E-06
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityY()
        {
            return new double[] {
                -5.30026E-08,
                -3.23281E-07,
                -1.00454E-06,
                -2.11387E-06,
                -3.49449E-06,
                -5.04341E-06,
                -6.74082E-06,
                -8.57872E-06,
                -1.05588E-05,
                -1.26869E-05
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.09, 0.09, 0.09
        private static double[] expectedSolidVelocityZ()
        {
            return new double[] {
                -3.31285E-05,
                -0.000132685,
                -0.000265916,
                -0.000399804,
                -0.000534109,
                -0.000668626,
                -0.000803237,
                -0.000937851,
                -0.001072391,
                -0.001206793,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //Great error values in comparison with comsol for the first 10 time steps. Results improve afterwards
        //GP : 0.05, 0.05, 0.05
        public static double[] expectedPressurevalues()
        {
            return new double[] {
                -2.62803E-05,
                -7.73618E-05,
                -0.000142106,
                -0.000178119,
                -0.000222339,
                -0.000265597,
                -0.000295448,
                -0.000282738,
                -0.00025486 ,
                -0.000182668
            };
        }
        
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dutdx_values()
        {
            return new double[] {
                3.09537E-06,
                1.86456E-05,
                5.75206E-05,
                0.000120403,
                0.000198256,
                0.000286194,
                0.000384862,
                0.000495531,
                0.000619937,
                0.000759913,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dvtdy_values()
        {
            return new double[] {
                6.36031E-06,
                3.87937E-05,
                0.000120545,
                0.000253664,
                0.000419339,
                0.000605209,
                0.000808898,
                0.001029447,
                0.001267052,
                0.001522424,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.09, 0.09, 0.09
        public static double[] expected_dwtdz_values()
        {
            return new double[] {
                -0.006214244,
                -0.024834619,
                -0.049594341,
                -0.07424962 ,
                -0.098808016,
                -0.123276063,
                -0.147651799,
                -0.17192928 ,
                -0.196099345,
                -0.220151411,
            };
        }

        public static double[] expected_div_vs_values()
        {
            return new double[10];
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdx_values()
        {
            return new double[] {
                -0.000510087,
                -0.000582082,
                0.000128969 ,
                0.005425463 ,
                0.011910756 ,
                0.020264921 ,
                0.028152223 ,
                0.033841904 ,
                0.038328132 ,
                0.039359541 ,
            };
        }

        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdy_values()
        {
            return new double[]
            {
                -0.001559998,
                -0.003938722,
                -0.006777791,
                -0.007825571,
                -0.010430786,
                -0.017679464,
                -0.026326288,
                -0.039003849,
                -0.052556255,
                -0.068290536,
            };
        }
        //MSOLVE reference values for Fz = {- 1e-4 / 4}
        //GP : 0.05, 0.05, 0.05
        public static double[] expected_dpdz_values()
        {
            return new double[]
            {
                -0.00085234,
                -0.001205138,
                -0.001344107,
                1.11914E-05,
                -0.000479235,
                -0.007219547,
                -0.016283798,
                -0.031968813,
                -0.049267171,
                -0.07062152,
            };
        }







    }
}
