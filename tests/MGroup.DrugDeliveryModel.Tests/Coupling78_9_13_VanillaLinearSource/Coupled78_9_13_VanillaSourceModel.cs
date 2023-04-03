using System;
using System.Collections.Generic;
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
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class Coupled78_9_13_VanillaSourceModel
    {
        public Eq78ModelProviderForStaggeredSolutionex7ref Eq78ModelProvider { get; set; }
        public CoxVanillaSourceModelBuilder CoxModelProvider { get; set; }
        public Eq9ModelProviderForStaggeredSolutionEx7Ref Eq9ModelProvider { get; set; }

        public Model[] model;



        //TODO put analysis time stepping where it belongs1 (pithanws sto Coupled7and9eqsSolution.cs h Coupled7and9eqsModel.cs)
        private GenericAnalyzerState[] analyzerStates, nlAnalyzerStates;
        private IParentAnalyzer[] parentAnalyzers;
        private IChildAnalyzer[] nlAnalyzers;
        private ISolver[] parentSolvers;


        public int CurrentTimeStep { get; set; }

        public GenericAnalyzerState[] AnalyzerStates => analyzerStates;
        public GenericAnalyzerState[] NLAnalyzerStates => nlAnalyzerStates;
        public IParentAnalyzer[] ParentAnalyzers => parentAnalyzers;
        public IChildAnalyzer[] NLAnalyzers => nlAnalyzers;
        public ISolver[] ParentSolvers => parentSolvers;
        public ComsolMeshReader Reader => reader;

        private ComsolMeshReader reader;

        private Dictionary<int, double> lambda;
        private Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints;
        private Dictionary<int, double[]> div_vs;
        private Dictionary<int, double[]> FluidSpeed;
        private double kth;

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public Coupled78_9_13_VanillaSourceModel(Eq78ModelProviderForStaggeredSolutionex7ref eq78ModelProvider, CoxVanillaSourceModelBuilder coxModelProvider,
                                     Eq9ModelProviderForStaggeredSolutionEx7Ref eq9ModelProvider, ComsolMeshReader comsolReader,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, Dictionary<int, double[]> FluidSpeed, double kth, double timeStep, double totalTime, int incrementsPerStep)
        {
            Eq9ModelProvider = eq9ModelProvider;
            CoxModelProvider = coxModelProvider;
            Eq78ModelProvider = eq78ModelProvider;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[3];
            nlAnalyzerStates = new GenericAnalyzerState[3];
            parentAnalyzers = new IParentAnalyzer[3];
            nlAnalyzers = new IChildAnalyzer[3];
            parentSolvers = new ISolver[3];

            reader = comsolReader;

            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            this.div_vs = div_vs;
            this.FluidSpeed = FluidSpeed;
            this.kth = kth;

            this.timeStep = timeStep;
            this.totalTime  = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[3];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            //---------------------------------------
            // WARNING: do not initialize shared dictionarys because they have been passed by refernce in ewuationModel bilders.
            //---------------------------------------


            // update Shared quantities of Coupled model
            //foreach (var elem in reader.ElementConnectivity)
            //{ 
            //    lambda[elem.Key]= lambda0;
            //}
            foreach (var elem in reader.ElementConnectivity)
            {
                pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
            }
            foreach (var elem in reader.ElementConnectivity)
            {
                div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
            }

            foreach (var elem in reader.ElementConnectivity)
            {//CALCULATE vf = kP + vs
                FluidSpeed[elem.Key] = new double[3]
                {
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][0] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][0]) * 10,
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][1] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][1]) * 10,
                    (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][2] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][2]) * 10,
                };
            }
            model = new Model[3];

            //Create model for eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddBoundaryConditions(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyper-elastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq13 (cox)
            model[2] = CoxModelProvider.GetModel();
            CoxModelProvider.AddBoundaryConditions(model[2]);
            (analyzers[2], solvers[2], nlAnalyzers[2]) = CoxModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }

        public void CreateModelFirstTime(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            if (!(CurrentTimeStep == 0))
            {
                foreach (var elem in reader.ElementConnectivity)
                {
                    pressureTensorDivergenceAtElementGaussPoints[elem.Key] = ((ConvectionDiffusionElement3D)model[0].ElementsDictionary[elem.Key]).pressureTensorDivergenceAtGaussPoints;
                }
                foreach (var elem in reader.ElementConnectivity)
                {
                    div_vs[elem.Key] = ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocityDivergence;
                }
                foreach (var elem in reader.ElementConnectivity)
                {//CALCULATE vf = kP + vs
                 //vs
                     FluidSpeed[elem.Key] = new double[3]
                     {
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][0] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][0]) * 10,
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][1] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][1]) * 10,
                         (-(pressureTensorDivergenceAtElementGaussPoints[elem.Key][0][2] * kth) + ((ContinuumElement3DGrowth)model[1].ElementsDictionary[elem.Key]).velocity[0][2]) * 10,
                     };
                }
            }


            model = new Model[3];
            
            //Create Initial Model eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddBoundaryConditions(model[0]);
            if(CurrentTimeStep==0)
            {
                //Eq78ModelProvider.AddEq78ModelInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyperelastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq8 (Cox)
            model[2] = CoxModelProvider.GetModel();
            CoxModelProvider.AddBoundaryConditions(model[2]);
            (analyzers[2], solvers[2], nlAnalyzers[2]) = CoxModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            for (int i = 0; i < analyzers.Length; i++)
            {
                analyzers[i].Initialize(true);
                if (analyzerStates[i] != null)
                {
                    analyzers[i].CurrentState = analyzerStates[i];
                }

                if (nlAnalyzerStates[i] != null)
                {
                    nlAnalyzers[i].CurrentState = nlAnalyzerStates[i];
                }
            }
        }

        public void SaveStateFromElements()
        {
            Eq9ModelProvider.SaveStateFromElements(model[1]);
        }
    }
}
