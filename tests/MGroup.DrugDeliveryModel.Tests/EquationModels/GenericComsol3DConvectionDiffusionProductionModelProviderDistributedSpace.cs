using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.MSolve.Discretization.Entities;
using MGroup.MSolve.Discretization;
using MGroup.FEM.Helpers;
using MGroup.DrugDeliveryModel.Tests.Commons;
using MGroup.Constitutive.Structural.InitialConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;

namespace MGroup.DrugDeliveryModel.Tests.EquationModels
{
    public class GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace
    {
        public Dictionary<int, double[]> ConvectionCoeffs;  //=> new[]  {1d, 1d, 1d};
        public double DiffusionCoeff;
        public Dictionary<int, double> DiffussionCoeffs;
        public Dictionary<int, double> DependentProductionCoeffs;
        public Dictionary<int, double> IndependentProductionCoeffs;
        private double CapacityCoeff;

        public ComsolMeshReader reader { get; private set; }

        public GenericComsol3DConvectionDiffusionProductionModelProviderDistributedSpace(ComsolMeshReader reader)
        {
            
            this.reader = reader;
        }

        public Model CreateModelFromComsolFile(Dictionary<int, double[]> convectionCoeffs,
            double diffusionCoeff, Dictionary<int, double> dependentProductionCoeffs,
            Dictionary<int, double> independentProductionCoeffs, double capacityCoeff,
            Func<double, double> productionFunc = null, Func<double, double> productionDeriv = null)
        {
            ConvectionCoeffs = convectionCoeffs;
            DiffusionCoeff = diffusionCoeff;
            DependentProductionCoeffs = dependentProductionCoeffs;
            IndependentProductionCoeffs = independentProductionCoeffs;
            CapacityCoeff = capacityCoeff;


            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);


            foreach (var node in reader.NodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }



            foreach (var elementConnectivity in reader.ElementConnectivity)
            {

                var material = new ConvectionDiffusionProperties(
                capacityCoeff: CapacityCoeff,
                diffusionCoeff: DiffusionCoeff,
                convectionCoeff: ConvectionCoeffs[elementConnectivity.Key],
                dependentSourceCoeff: productionDeriv != null ? 0 : DependentProductionCoeffs[elementConnectivity.Key],
                independentSourceCoeff: IndependentProductionCoeffs[elementConnectivity.Key]);

                var elementFactory = new ConvectionDiffusionElement3DFactory(material);
                var element = elementFactory.CreateElement(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2);
                if (productionDeriv != null)
                {
                    element.LinearProduction = false;
                    element.ProductionFunction = productionFunc;
                    element.ProductionFunctionDerivative = productionDeriv;

                }
                element.ID = elementConnectivity.Key; ;
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }


            return model;
        }
        public Model CreateModelFromComsolFile(Dictionary<int, double[]> convectionCoeffs,
            double diffusionCoeff, Dictionary<int, double> dependentProductionCoeffs,
            Dictionary<int, double> independentProductionCoeffs, double capacityCoeff,
             Dictionary<int, Func<double, double>> productionFunc, Dictionary<int, Func<double, double>> productionDeriv)
        {
            ConvectionCoeffs = convectionCoeffs;
            DiffusionCoeff = diffusionCoeff;
            DependentProductionCoeffs = dependentProductionCoeffs;
            IndependentProductionCoeffs = independentProductionCoeffs;
            CapacityCoeff = capacityCoeff;


            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);


            foreach (var node in reader.NodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }



            foreach (var elementConnectivity in reader.ElementConnectivity)
            {

                var material = new ConvectionDiffusionProperties(
                capacityCoeff: CapacityCoeff,
                diffusionCoeff: DiffusionCoeff,
                convectionCoeff: ConvectionCoeffs[elementConnectivity.Key],
                dependentSourceCoeff: productionDeriv != null ? 0 : DependentProductionCoeffs[elementConnectivity.Key],
                independentSourceCoeff: IndependentProductionCoeffs[elementConnectivity.Key]);

                var elementFactory = new ConvectionDiffusionElement3DFactory(material);
                var element = elementFactory.CreateElement(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2);
                if (productionDeriv != null)
                {
                    element.LinearProduction = false;
                    element.ProductionFunction = productionFunc[elementConnectivity.Key];
                    element.ProductionFunctionDerivative = productionDeriv[elementConnectivity.Key];

                }
                element.ID = elementConnectivity.Key; ;
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }


            return model;
        }

        public Model CreateModelFromComsolFile(Dictionary<int, double[]> convectionCoeffs,
            Dictionary<int, double> difussionCoeffs, Dictionary<int, double> dependentProductionCoeffs,
            Dictionary<int, double> independentProductionCoeffs, double capacityCoeff,
            Func<double, double> productionFunc = null, Func<double, double> productionDeriv = null)
        {
            ConvectionCoeffs = convectionCoeffs;
            DiffussionCoeffs = difussionCoeffs;
            DependentProductionCoeffs = dependentProductionCoeffs;
            IndependentProductionCoeffs = independentProductionCoeffs;
            CapacityCoeff = capacityCoeff;


            var model = new Model();
            model.SubdomainsDictionary[0] = new Subdomain(id: 0);


            foreach (var node in reader.NodesDictionary.Values)
            {
                model.NodesDictionary.Add(node.ID, node);
            }



            foreach (var elementConnectivity in reader.ElementConnectivity)
            {

                var material = new ConvectionDiffusionProperties(
                capacityCoeff: CapacityCoeff,
                diffusionCoeff: DiffussionCoeffs[elementConnectivity.Key],
                convectionCoeff: ConvectionCoeffs[elementConnectivity.Key],
                dependentSourceCoeff: productionDeriv != null ? 0 : DependentProductionCoeffs[elementConnectivity.Key],
                independentSourceCoeff: IndependentProductionCoeffs[elementConnectivity.Key]);

                var elementFactory = new ConvectionDiffusionElement3DFactory(material);
                var element = elementFactory.CreateElement(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2);
                if (productionDeriv != null)
                {
                    element.LinearProduction = false;
                    element.ProductionFunction = productionFunc;
                    element.ProductionFunctionDerivative = productionDeriv;

                }
                element.ID = elementConnectivity.Key; ;
                model.ElementsDictionary.Add(elementConnectivity.Key, element);
                model.SubdomainsDictionary[0].Elements.Add(element);
            }


            return model;
        }


        /// <summary>
        /// gia to model workingtetmesh155 to 0.1, 0.1, 0.,1 kuvo 
        /// </summary>
        /// <param name="model"></param>
        /// <param name="maxZ">value =2</param>
        /// <param name="topValueprescribed"></param>
        /// <param name="minZ">value = 0</param>
        /// <param name="bottomValueprescribed"></param>
        /// <returns></returns>
        public void AddTopAndBottomBCs(Model model, double modelMaxZ, double topValueprescribed, double modelMinZ, double bottomValueprescribed)
        {
            
            var topNodes = new List<INode>();
            var bottomNodes = new List<INode>();
            var innerBulkNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxZ - node.Z) < 1E-9) topNodes.Add(node);
                else if (Math.Abs(modelMinZ - node.Z) < 1E-9) bottomNodes.Add(node);
                else innerBulkNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in topNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, topValueprescribed));            
            }
            foreach (var node in bottomNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, bottomValueprescribed));
            }

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }

        public void AddZplusBCs(Model model, double modelMaxZ, double topValueprescribed)
        {

            var topNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxZ - node.Z) < 1E-9) topNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in topNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, topValueprescribed));
            }
            

            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }

        public void AddXplusBCs(Model model, double modelMaxX, double ValueprescribedformaxXcoord)
        {

            var maxXCoordNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxX - node.X) < 1E-9) maxXCoordNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in maxXCoordNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, ValueprescribedformaxXcoord));
            }


            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }

        public void AddYplusBCs(Model model, double modelMaxX, double ValueprescribedformaxXcoord)
        {

            var maxXCoordNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxX - node.X) < 1E-9) maxXCoordNodes.Add(node);
            }

            var dirichletBCs = new List<NodalUnknownVariable>();
            foreach (var node in maxXCoordNodes)
            {
                dirichletBCs.Add(new NodalUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, ValueprescribedformaxXcoord));
            }


            model.BoundaryConditions.Add(new ConvectionDiffusionBoundaryConditionSet(
                dirichletBCs,
                new INodalConvectionDiffusionNeumannBoundaryCondition[] { }
            ));

        }



        public void AddInitialConditionsForTheRestOfBulkNodesMInusTopAndBottom(Model model, double modelMaxZ, double modelMinZ, double prescribedInitialConditionValueForBulk)
        {

            var topNodes = new List<INode>();
            var bottomNodes = new List<INode>();
            var innerBulkNodes = new List<INode>();
            foreach (var node in model.NodesDictionary.Values)
            {
                if (Math.Abs(modelMaxZ - node.Z) < 1E-9) topNodes.Add(node);
                else if (Math.Abs(modelMinZ - node.Z) < 1E-9) bottomNodes.Add(node);
                else innerBulkNodes.Add(node);
            }

            var intitalConditions = new List<INodalConvectionDiffusionInitialCondition>();
            foreach (var node in innerBulkNodes)
            {
                if ((!topNodes.Contains(node)) && (!bottomNodes.Contains(node)))
                {
                    intitalConditions.Add(new NodalInitialUnknownVariable(node, ConvectionDiffusionDof.UnknownVariable, prescribedInitialConditionValueForBulk));
                }

            }


            model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(intitalConditions,
                new DomainInitialUnknownVariable[]
                { }));

            
        }



    }
}
