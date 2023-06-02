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
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Solution;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using System.Xml.Linq;
using TriangleNet.Tools;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class CopledUnewmarkPnewmarkT
    {
        public Eq78ModelProviderForStaggeredSolutionex7ref eq78ModelProviderForCouplin { get; }
        public Eq9ModelProvider eq9ModelProvider { get; }
        public TCellModelProvider TCellModelProvider { get; }

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
        private Dictionary<int, double[][]> SolidVelocityAtElementGaussPoints;

        private double timeStep;
        private double totalTime;

        private int incrementsPerStep;

        public CopledUnewmarkPnewmarkT(Eq78ModelProviderForStaggeredSolutionex7ref Eq78ModelProviderForStaggeredSolutionex7ref,
                                                 Eq9ModelProvider solidPhaseProvider,
                                                 TCellModelProvider tCellModelProvider,
                                                 ComsolMeshReader comsolReader,
            Dictionary<int, double> lambda, Dictionary<int, double[][]> pressureTensorDivergenceAtElementGaussPoints,
            Dictionary<int, double[]> div_vs, Dictionary<int, double[][]> velocityAtElementGaussPoints, double timeStep, double totalTime, int incrementsPerStep)
        {
            eq9ModelProvider = solidPhaseProvider;
            eq78ModelProviderForCouplin = Eq78ModelProviderForStaggeredSolutionex7ref;
            TCellModelProvider = tCellModelProvider;
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
            
            this.SolidVelocityAtElementGaussPoints = velocityAtElementGaussPoints;
            
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
          

            eq78ModelProviderForCouplin.UpdatePressureDivergenceDictionary(pressureTensorDivergenceAtElementGaussPoints, ParentSolvers[0], NLAnalyzers[0], model[0], eq78ModelProviderForCouplin.algebraicModel);

            (ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();
            
            model[1].BoundaryConditions.Clear();
            var velocities = eq9ModelProvider.GetVelocities2();

            foreach (var elem in reader.ElementConnectivity)
            {
                SolidVelocityAtElementGaussPoints[elem.Key] = new double[][]{new double[]{
                    velocities[elem.Key][0][0] * 1000,
                    velocities[elem.Key][0][1] * 1000,
                    velocities[elem.Key][0][2] * 1000 }};
            }

            var velocityDIVs = eq9ModelProvider.GetVelocityDIV();
            div_vs = velocityDIVs;
            eq78ModelProviderForCouplin.div_vs = div_vs;
            #region analyzers are cleared here as they are overwritten in next commands
            #endregion

            model = new Model[3];
            
            //Create model for eq78 (fluid pressure)
            model[0] = eq78ModelProviderForCouplin.GetModel();
            eq78ModelProviderForCouplin.AddBoundaryConditions(model[0]);
            (analyzers[0], solvers[0], nlAnalyzers[0]) = eq78ModelProviderForCouplin.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            //Create model for eq9 (hyper-elastic material)
            model[1] = eq9ModelProvider.GetModel();
            eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);
            ParentAnalyzers[1] = analyzers[1];


            model[2] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[2]);
            (analyzers[2], solvers[2], nlAnalyzers[2]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep);

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
                eq78ModelProviderForCouplin.UpdatePressureDivergenceDictionary(pressureTensorDivergenceAtElementGaussPoints, ParentSolvers[0], NLAnalyzers[0], model[0], eq78ModelProviderForCouplin.algebraicModel);

                model[1].BoundaryConditions.Clear();
                var velocities = eq9ModelProvider.GetVelocities2();

                foreach (var elem in reader.ElementConnectivity)
                {
                    SolidVelocityAtElementGaussPoints[elem.Key] = new double[][]{new double[]{
                    velocities[elem.Key][0][0] * 1000,
                    velocities[elem.Key][0][1] * 1000,
                    velocities[elem.Key][0][2] * 1000 }};
                }

                var velocityDIVs = eq9ModelProvider.GetVelocityDIV();
                div_vs = velocityDIVs;
                eq78ModelProviderForCouplin.div_vs = div_vs;
                #region analyzers are cleared here as they are overwritten in next commands
                #endregion

            }


            model = new Model[3];
            
            model[0] = eq78ModelProviderForCouplin.GetModel();
            eq78ModelProviderForCouplin.AddBoundaryConditions(model[0]);
            if(CurrentTimeStep==0)
            {
                //Eq78ModelProvider.AddEq78ModelInitialConditions(model[0]);
            }
            (analyzers[0], solvers[0], nlAnalyzers[0]) = eq78ModelProviderForCouplin.GetAppropriateSolverAnalyzerAndLog(model[0], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);

            model[1] = eq9ModelProvider.GetModel();
            eq9ModelProvider.AddBoundaryConditions(model[1]);
            (analyzers[1], solvers[1], nlAnalyzers[1]) = eq9ModelProvider.GetAppropriateSolverAnalyzerAndLog(model[1], timeStep, totalTime, CurrentTimeStep, incrementsPerStep);
            ParentAnalyzers[1] = analyzers[1];

            model[2] = TCellModelProvider.GetModel();
            TCellModelProvider.AddBoundaryConditions(model[2]);
            if (CurrentTimeStep == 0)
            {
                TCellModelProvider.AddInitialConditions(model[2]);
            }
            (analyzers[2], solvers[2], nlAnalyzers[2]) = TCellModelProvider.GetAppropriateSolverAnalyzerAndLog(model[2], timeStep, totalTime, CurrentTimeStep);

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
            eq9ModelProvider.SaveStateFromElements(model[1]);
        }

        //private void UpdatePressureAndGradients(double[] localDisplacements)
        //{
        //    IReadOnlyList<Matrix> shapeFunctionNaturalDerivatives;
        //    shapeFunctionNaturalDerivatives = Interpolation.EvaluateNaturalGradientsAtGaussPoints(QuadratureForConsistentMass);
        //    var jacobians = shapeFunctionNaturalDerivatives.Select(x => new IsoparametricJacobian3D(Nodes, x));
        //    Matrix[] jacobianInverse = jacobians.Select(x => x.InverseMatrix.Transpose()).ToArray();
        //    for (int gp = 0; gp < QuadratureForConsistentMass.IntegrationPoints.Count; ++gp)
        //    {
        //        double[] dphi_dnatural = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
        //        for (int i1 = 0; i1 < shapeFunctionNaturalDerivatives[gp].NumRows; i1++)
        //        {
        //            dphi_dnatural[0] += shapeFunctionNaturalDerivatives[gp][i1, 0] * localDisplacements[i1];
        //            dphi_dnatural[1] += shapeFunctionNaturalDerivatives[gp][i1, 1] * localDisplacements[i1];
        //            dphi_dnatural[2] += shapeFunctionNaturalDerivatives[gp][i1, 2] * localDisplacements[i1];
        //        }

        //        var dphi_dnaturalMAT = Matrix.CreateFromArray(dphi_dnatural, 1, 3);

        //        var dphi_dcartesian = dphi_dnaturalMAT * jacobianInverse[gp].Transpose();

        //        pressureTensorDivergenceAtGaussPoints[gp] = new double[3] { dphi_dcartesian[0, 0], dphi_dcartesian[0, 1], dphi_dcartesian[0, 2] };




        //    }
        //}
    }
}