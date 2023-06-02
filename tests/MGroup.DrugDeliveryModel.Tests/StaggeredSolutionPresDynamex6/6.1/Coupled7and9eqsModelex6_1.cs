using System.Collections.Generic;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Interpolation.Jacobians;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Numerics.Integration;
using MGroup.NumericalAnalyzers.Dynamic;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled7and9eqsModelex6_1
    {
        public Eq78ModelProviderForStaggeredSolutionex6_1 Eq78ModelProvider { get; set; }
        public Eq9ModelProviderForStaggeredSolutionex6_1 Eq9ModelProvider { get; set; }

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

        private Dictionary<int, double[][]> SolidVelocityAtElementGaussPoints = new Dictionary<int, double[][]>();
        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        int nGaussPoints = 1;
        static double initial_dp_dx = 0.0;
        static double initial_dp_dy = 0.0;
        static double initial_dp_dz = 0.0;
        static double velocityDivInitialVal = 0;
        public Coupled7and9eqsModelex6_1(Eq78ModelProviderForStaggeredSolutionex6_1 eq78ModelProvider,
                                     Eq9ModelProviderForStaggeredSolutionex6_1 eq9ModelProvider, ComsolMeshReader comsolReader,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, double timeStep, double totalTime, int incrementsPerStep)
        {
            Eq9ModelProvider = eq9ModelProvider;
            Eq78ModelProvider = eq78ModelProvider;
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;


            analyzerStates = new GenericAnalyzerState[2];
            nlAnalyzerStates = new GenericAnalyzerState[2];
            parentAnalyzers = new IParentAnalyzer[2];
            nlAnalyzers = new IChildAnalyzer[2];
            parentSolvers = new ISolver[2];

            reader = comsolReader;

            this.pressureTensorDivergenceAtElementGaussPoints = pressureTensorDivergenceAtElementGaussPoints;
            this.lambda = lambda;
            this.div_vs = div_vs;

            this.timeStep = timeStep;
            this.totalTime = totalTime;
            this.incrementsPerStep = incrementsPerStep;

            // intialize array ofm models1.
            model = new Model[2];
        }


        public void CreateModel(IParentAnalyzer[] analyzers, ISolver[] solvers)
        {
            Eq78ModelProvider.UpdatePressureDivergenceDictionary(pressureTensorDivergenceAtElementGaussPoints, ParentSolvers[0], NLAnalyzers[0], model[0], Eq78ModelProvider.algebraicModel);
            (ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();
            model[1].BoundaryConditions.Clear();
            var velocities = Eq9ModelProvider.GetVelocities2();

            foreach (var elem in reader.ElementConnectivity)
            {
                SolidVelocityAtElementGaussPoints[elem.Key] = new double[][]{new double[]{
                    velocities[elem.Key][0][0] * 1000,
                    velocities[elem.Key][0][1] * 1000,
                    velocities[elem.Key][0][2] * 1000 }};
            }

            var velocityDIVs = Eq9ModelProvider.GetVelocityDIV();
            div_vs = velocityDIVs;
            Eq78ModelProvider.div_vs = div_vs;

            model = new Model[2];

            //Create model for eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddAllBoundaryNodesBC(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyper-elastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBottomLeftFrontBackBCs(model[1]);
            Eq9ModelProvider.AddEq9ModelLoads(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

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
                Eq78ModelProvider.UpdatePressureDivergenceDictionary(pressureTensorDivergenceAtElementGaussPoints, ParentSolvers[0], NLAnalyzers[0], model[0], Eq78ModelProvider.algebraicModel);

                model[1].BoundaryConditions.Clear();
                var velocities = Eq9ModelProvider.GetVelocities2();

                foreach (var elem in reader.ElementConnectivity)
                {
                    SolidVelocityAtElementGaussPoints[elem.Key] = new double[][]{new double[]{
                    velocities[elem.Key][0][0] * 1000,
                    velocities[elem.Key][0][1] * 1000,
                    velocities[elem.Key][0][2] * 1000 }};
                }

                var velocityDIVs = Eq9ModelProvider.GetVelocityDIV();
                div_vs = velocityDIVs;
                Eq78ModelProvider.div_vs = div_vs;
            }
            model = new Model[2];

            //Create Initial Model eq78 (fluid pressure)
            model[0] = Eq78ModelProvider.GetModel();
            Eq78ModelProvider.AddAllBoundaryNodesBC(model[0]);
            if (CurrentTimeStep == 0)
            {
                //Eq78ModelProvider.AddEq78ModelInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = Eq78ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyperelastic material)
            model[1] = Eq9ModelProvider.GetModel();
            Eq9ModelProvider.AddBottomLeftFrontBackBCs(model[1]);
            Eq9ModelProvider.AddEq9ModelLoads(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = Eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

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




    }
}
