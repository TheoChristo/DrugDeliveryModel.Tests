using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.DrugDeliveryModel.Tests.EquationModels;
using MGroup.DrugDeliveryModel.Tests.Integration;
using MGroup.DrugDeliveryModel.Tests.Materials;
using MGroup.DrugDeliveryModel.Tests.PreliminaryModels;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.AnalysisWorkflow;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Numerics.Integration;
using MGroup.MSolve.Solution.AlgebraicModel;
using MGroup.NumericalAnalyzers;
using MGroup.NumericalAnalyzers.Discretization.NonLinear;
using MGroup.NumericalAnalyzers.Dynamic;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.NumericalAnalyzers.Staggered;
using MGroup.Solvers.Direct;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection.PortableExecutable;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using Xunit;

namespace MGroup.DrugDeliveryModel.Tests.IntegrationTests1
{
    public class HexaForcesTest
    {
        private const double timeStep = 0.00001; // in sec
        const double totalTime = 0.01; // in sec
        static int incrementsPertimeStep = 1;
        static int currentTimeStep = 0;

        static double density = 1;

        static double miTumor = 22.44; //Kpa
        static double kappaTumor = 216.7; //Kpa

        static double lambda0 = 1;
        [Theory]

        [InlineData("../../../DataFiles/simple2HexaMesh.mphtxt")]
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

            var constrainedNodeIds = new[] {
                Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new[] {0d,0d,0d}, 1e-2),
                Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new[] {0d,1d,0d}, 1e-2),
                Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new[] {1d,0d,0d}, 1e-2),
                Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, new[] {1d,1d,0d}, 1e-2),
            };


            var nodes = comsolReader.NodesDictionary;
            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);

            foreach (var node in nodes.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }

            var materialTumor = new NeoHookeanMaterial3dJ3Isochoric(miTumor, kappaTumor);

            var elasticMaterial = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
            var dynamicMaterial = new TransientAnalysisProperties(density: density, rayleighCoeffMass: 0, rayleighCoeffStiffness: 0);
            var elementFactory = new ContinuumElement3DFactory(elasticMaterial, dynamicMaterial);

            var dirichletBCs = new INodalDisplacementBoundaryCondition[constrainedNodeIds.Length*3];
            for (int i = 0; i < constrainedNodeIds.Length; i++)
            {
                dirichletBCs[3*i] = new NodalDisplacement(model.NodesDictionary[constrainedNodeIds[i]], StructuralDof.TranslationX, 0);
                dirichletBCs[3*i+1] = new NodalDisplacement(model.NodesDictionary[constrainedNodeIds[i]], StructuralDof.TranslationY, 0);
                dirichletBCs[3*i+2] = new NodalDisplacement(model.NodesDictionary[constrainedNodeIds[i]], StructuralDof.TranslationZ, 0);
            }


            var volumeLoads = new List<IElementStructuralNeumannBoundaryCondition>();

            //var domains = new Dictionary<int, double[]>(2);
            foreach (var elementConnectivity in comsolReader.ElementConnectivity)
            {
                var element = elementFactory.CreateNonLinearElementGrowt(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2, materialTumor, dynamicMaterial, lambda[elementConnectivity.Key]);
                //element.volumeForce = pressureTensorDivergenceAtElementGaussPoints[elementConnectivity.Key][0];

                Console.WriteLine("Volume force on element:" + element.ID);
                var volumeForceX = new ElementDistributedLoad(element, StructuralDof.TranslationX, -2);
                var volumeForceY = new ElementDistributedLoad(element, StructuralDof.TranslationY, -9);
                var volumeForceZ = new ElementDistributedLoad(element, StructuralDof.TranslationZ, -5);
                volumeLoads.AddRange(new[] { volumeForceX, volumeForceY, volumeForceZ });

                element.ID = elementConnectivity.Key;
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }

            var modelNeumannConditions = new StructuralBoundaryConditionSet(dirichletBCs, null, null, volumeLoads, null, null);
            model.BoundaryConditions.Add(modelNeumannConditions);

            var solverFactory = new SkylineSolver.Factory() { FactorizationPivotTolerance = 1e-8 };
            var algebraicModel = solverFactory.BuildAlgebraicModel(model);
            var solver = solverFactory.BuildSolver(algebraicModel);
            var provider = new ProblemStructural(model, algebraicModel);
            var loadControlAnalyzerBuilder = new LoadControlAnalyzer.Builder(algebraicModel, solver, provider, numIncrements: incrementsPertimeStep)
            {
                ResidualTolerance = 1E-8,
                MaxIterationsPerIncrement = 100,
                NumIterationsForMatrixRebuild = 1
            };
            var nlAnalyzer = loadControlAnalyzerBuilder.Build();
            var loadControlAnalyzer = (LoadControlAnalyzer)nlAnalyzer;

            var staticAnalyzer = new StaticAnalyzer(algebraicModel, provider, loadControlAnalyzer);

            
            staticAnalyzer.Initialize();
            var rhs = provider.GetRhs();
            var rhsExpected = expectedRhs();
            //staticAnalyzer.Solve();

            //var analyzerBuilder = new NewmarkDynamicAnalyzer.Builder(algebraicModel, provider, loadControlAnalyzer, timeStep: timeStep, totalTime: totalTime, false, currentStep: currentTimeStep);
            //analyzerBuilder.SetNewmarkParametersForConstantAcceleration();
            //var analyzer = analyzerBuilder.Build();
            //analyzer.Initialize();
            //analyzer.Solve();
        }

        public static double[] expectedRhs()
        {
            return new double[] {
             0.00,
             0.00,
             -1.25,
             0.00,
             0.00,
             -1.25,
             0.00,
             0.00,
             -1.25,
             0.00,
             0.00,
             -0.625,
             0.00,
             0.00,
             -1.25,
             0.00,
             0.00,
             -0.625,
             0.00,
             0.00,
             -0.625,
             0.00,
             0.00,
             -0.625 };
        }

    }
}
