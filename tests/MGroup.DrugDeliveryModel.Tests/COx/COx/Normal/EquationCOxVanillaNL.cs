using BC = MGroup.DrugDeliveryModel.Tests.Commons.BoundaryAndInitialConditionsUtility.BoundaryConditionCase;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
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
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using System.Diagnostics;
using MGroup.FEM.ConvectionDiffusion.Tests.Commons;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class EquationCOxVanillaNL
    {
        //---------------------------------------Equation Cox---------------------------------------

        //---------------------------------------Variables------------------------------------------

        private Dictionary<int, double[]> DomainFluidVelocity = new Dictionary<int, double[]>(); // 2.32E-4 [m/s]

        //private Dictionary<int, double> DomainIndependentSource = new Dictionary<int, double>();

        //private Dictionary<int, double> DomainDependentSource = new Dictionary<int, double>(); 

        double FluidInit = -2.32;
        //double FluidInit = -1.0;
        /// <summary>
        /// Diffusivity of oxygen [m2/s]
        /// </summary>
        private const double Dox = 1.79E-4; // [m2/s]

        /// <summary>
        /// Oxygen uptake [mol/(m3*s)]
        /// </summary>
        private const double Aox = 2.5463E-2; // [mol/(m3*s)]

        /// <summary>
        /// Oxygen uptake [mol/m3]
        /// </summary>
        private const double Kox = 4.64E-3; // [mol / m3]

        /// <summary>
        /// Oxygen permeability across tumor vessel walls [m/s]
        /// </summary>
        private const double PerOx = 3.55E-4; // [m/s]

        /// <summary>
        /// Vascular Density [1/m]
        /// </summary>
        private const double Sv = 7E3; // [1/m]

        /// <summary>
        /// Initial Oxygen Concentration [mol/m3]
        /// </summary>
        private const double CInitialOx = 0; // [mol/m3]

        /// <summary>
        /// Rhs [mol/m3]
        /// </summary>
        private const double CiOx = 0.2; // [mol/m3]
        const double TInit = 500;
        /// <summary>
        /// Cancer cell density [1]
        /// </summary>
        private Dictionary<int, double> T = new Dictionary<int, double>();// 500 [cells]

        //---------------------------------------Logging----------------------------------
        /// <summary>
        /// The degree of freedom that will be monitored for equation cox
        /// </summary>
        private ConvectionDiffusionDof coxMonitorDOF = ConvectionDiffusionDof.UnknownVariable;

        /// <summary>
        /// The coordinates of the monitored node
        /// </summary>
        private double[] monitorNodeCoords = { 0.09, 0.0, 0.08 };


        //---------------------------------------Time Discretization Specs------------------------------
        private const double TotalTime = 20E-5;

        /// <summary>
        /// For increased accuracy use time-step of order 1E-5
        /// </summary>
        private const double TimeStep = 1E-5;// sec

        /// <summary>
        /// Simplified version of the independent production term without non-linear term
        /// </summary>
        public Func<double> IndependentLinearSource = () => PerOx * Sv * CiOx;

        public Func<double> DependentLinearSource => null;// () => -PerOx * Sv;

        private Dictionary<int, Func<double, double>> ProductionFuncsWithoutConstantTerm = new Dictionary<int, Func<double, double>>();
        public Func<double, double> getProductionFuncWithoutConstantTerm(int i)
        {
            return (double Cox) => -PerOx * Sv * (CiOx - Cox) - Aox * T[i] * Cox / (Cox + Kox);
            //return (double Cox) => -PerOx * Sv * Cox; //Linear
        }

        private Dictionary<int, Func<double, double>> ProductionFuncsWithoutConstantTermDerivative = new Dictionary<int, Func<double, double>>();
        public Func<double, double> getProductionFuncWithoutConstantTermDerivative(int i)
        {
            return (double Cox) => -PerOx * Sv - Aox * T[i] / (Cox + Kox) + Aox * T[i] * Cox * Math.Pow(Cox + Kox, -2);
            //return (double Cox) => -PerOx * Sv; //Linear
        }

        static double[] expectedLinSolution = new double[]
            {0.00027617099601129385,
             0.00082795347593641352,
             0.0013786176765330069,
             0.0019281650868235782,
             0.0024765971944125868,
             0.00302391548549161,
             0.003570121444844468,
             0.0041152165558522892,
             0.0046592023004985586,
             0.0052020801593741341,
             0.0057438516116821775,
             0.0062845181352431112,
             0.0068240812064995068,
             0.007362542300520903,
             0.0078999028910086543,
             0.0084361644503007076,
             0.0089713284493762886,
             0.0095053963578606761,
             0.010038369644029817,
             0.01057024977481494};

        static double[] expectedFromComsol = new double[]
        {0.0007129565223671997, 0.0010666609276386, 0.0014182726753743, 0.0017679034690053, 0.0021155406016867, 0.0024612345221367, 0.0028050506063947, 0.0031470540006927, 0.0034873053358745, 0.0038258599357506, 0.0041627680836797, 0.0044980755555751, 0.00483182417161, 0.0051640522969102, 0.0054947952781904, 0.0058240858200981, 0.0061519543090138, 0.0064784290921262, 0.0068035367186388, 0.0071273021488608};

        public void EquationsTests13DistributedModelBuilder()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }
        
        [Theory]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]
        public void SolveEquationCOxNonLinearProduction(string fileName)
        {
            var mesh = new ComsolMeshReader(fileName);
            foreach (var element in mesh.ElementConnectivity)
            {
                DomainFluidVelocity.Add(element.Key, new double[] { FluidInit, FluidInit, FluidInit });
                if (DependentLinearSource == null)
                {
                    T.Add(element.Key, TInit);
                    ProductionFuncsWithoutConstantTerm.Add(element.Key, getProductionFuncWithoutConstantTerm(element.Key));
                    ProductionFuncsWithoutConstantTermDerivative.Add(element.Key, getProductionFuncWithoutConstantTermDerivative(element.Key));
                }
            }
            var convectionDiffusionDirichletBC =
            new List<(BC, ConvectionDiffusionDof[], double[][], double[])>()
            {
                (BC.TopRightBackDiriclet, new ConvectionDiffusionDof[] { ConvectionDiffusionDof.UnknownVariable }, new double[2][]{new double[3] {0,0,0},new double[3] {0.1,0.1,0.1}}, new double[]{0.2}),


            };
            var convectionDiffusionNeumannBC = new List<(BoundaryAndInitialConditionsUtility.BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])>();

            var nodeIdToMonitor = Utilities.FindNodeIdFromNodalCoordinates(mesh.NodesDictionary, monitorNodeCoords, 1e-4);

            var modelBuilder = new CoxModelBuilder(mesh, DomainFluidVelocity,
                                                    Dox, Aox, Kox, PerOx, Sv, CiOx, T, CInitialOx,
                                                    IndependentLinearSource, DependentLinearSource, ProductionFuncsWithoutConstantTerm, ProductionFuncsWithoutConstantTermDerivative,
                                                    nodeIdToMonitor, coxMonitorDOF,
                                                    convectionDiffusionDirichletBC, convectionDiffusionNeumannBC);
            var model = modelBuilder.GetModel();
            modelBuilder.AddBoundaryConditions(model);

            (var analyzer, var solver, var nlAnalyzers) = modelBuilder.GetAppropriateSolverAnalyzerAndLog(model, TimeStep, TotalTime, 0, 1);

            ((NewmarkDynamicAnalyzer)analyzer).ResultStorage = new ImplicitIntegrationAnalyzerLog();

            analyzer.Initialize(true);
            analyzer.Solve();

            int totalNewmarkstepsNum = (int)Math.Truncate(TotalTime / TimeStep);
            var cox = new double[totalNewmarkstepsNum];
            for (int i1 = 0; i1 < totalNewmarkstepsNum; i1++)
            {
                var timeStepResultsLog = ((NewmarkDynamicAnalyzer)analyzer).ResultStorage.Logs[i1];
                cox[i1] = ((DOFSLog)timeStepResultsLog).DOFValues[model.GetNode(nodeIdToMonitor), coxMonitorDOF];
            }
            
            CSVExporter.ExportVectorToCSV(cox, "../../../Integration/coxNL.csv");
            Assert.True(ResultChecker.CheckResults(cox, expectedFromComsol, 1E-6));

            //Console.WriteLine("FINISHED solving Cox Non-Linear prod");

        }



    }
}
