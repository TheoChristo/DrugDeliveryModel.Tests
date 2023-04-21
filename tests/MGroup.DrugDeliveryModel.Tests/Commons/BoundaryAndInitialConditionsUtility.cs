using System;
using System.Collections.Generic;
using System.Linq;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.ConvectionDiffusion.BoundaryConditions;
using MGroup.Constitutive.ConvectionDiffusion.InitialConditions;
using MGroup.Constitutive.Structural;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.MSolve.Discretization.BoundaryConditions;
using MGroup.MSolve.Discretization.Entities;

namespace MGroup.DrugDeliveryModel.Tests.Commons;

public class BoundaryAndInitialConditionsUtility
{
	public enum BoundaryConditionCase
	{
		LeftDirichlet,
		LeftPointFlux,
		LeftDistributedFlux,
		
		RightDirichlet,
		RightPointFlux,
		RightDistributedFlux,
		
		FrontDirichlet,
		FrontPointFlux,
		FrontDistributedFlux,
		
		BackDirichlet,
		BackPointFlux,
		BackDistributedFlux,
		
		BottomDirichlet,
		BottomPointFlux,
		BottomDistributedFlux,
		
		TopDirichlet,
		TopPointFlux,
		TopDistributedFlux,

		TopRightBackDiriclet
	}
	
	public static void AssignStructuralDirichletBCsToModel(Model model,
		List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
		double tolerance)
	{
		var constraints = new List<INodalDisplacementBoundaryCondition>();
		
		foreach (var bcSet in bcs)
		{
			var bcFace = bcSet.Item1;
			var bcDofs = bcSet.Item2;
			var faceCoords = bcSet.Item3;
			var bcValues = bcSet.Item4;
			
			foreach (var node in model.NodesDictionary)
			{
				var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };

				switch (bcFace)
				{
					case BoundaryConditionCase.RightDirichlet or BoundaryConditionCase.LeftDirichlet:
					{
						if (Math.Abs(nodalCoords[0] - faceCoords[0][0]) < tolerance)
						{
							AddBoundedDisplacementDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
					case BoundaryConditionCase.FrontDirichlet or BoundaryConditionCase.BackDirichlet:
					{
						if (Math.Abs(nodalCoords[1] - faceCoords[0][1]) < tolerance)
						{
							AddBoundedDisplacementDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
					case BoundaryConditionCase.TopDirichlet or BoundaryConditionCase.BottomDirichlet:
					{
						if (Math.Abs(nodalCoords[2] - faceCoords[0][2]) < tolerance)
						{
							AddBoundedDisplacementDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
				}
			}
		}

		var modelDirichletBCs =
			new StructuralBoundaryConditionSet(constraints, new List<INodalLoadBoundaryCondition>());
		model.BoundaryConditions.Add(modelDirichletBCs);
	}

	public static void AssignStructuralNeumannBCsToModel(Model model, List<(BoundaryConditionCase, StructuralDof[], double[][], double[])> bcs,
	    double tolerance)
    {
        var loads = new List<INodalLoadBoundaryCondition>();
        
        foreach (var node in model.NodesDictionary)
        {
        	var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };

            foreach (var bcSet in bcs)
            {
	            //TODO : Should not delete cases. Point loads is one and distributed loads is another.
                var loadCase = bcSet.Item1;
                var bcDofs = bcSet.Item2;
                var loadCoords = bcSet.Item3;
                var bcValues = bcSet.Item4;

                foreach (var loadCoordinates in loadCoords)
                {
	                switch (loadCase)
	                {
		                case BoundaryConditionCase.LeftPointFlux or BoundaryConditionCase.RightPointFlux or 
			                BoundaryConditionCase.FrontPointFlux or BoundaryConditionCase.BackPointFlux or 
			                BoundaryConditionCase.TopPointFlux or BoundaryConditionCase.BottomPointFlux:
		                {
			                if (Math.Abs(nodalCoords[0] - loadCoordinates[0]) < tolerance &&
			                    Math.Abs(nodalCoords[1] - loadCoordinates[1]) < tolerance &&
			                    Math.Abs(nodalCoords[2] - loadCoordinates[2]) < tolerance)
			                {
				                AddNodalForceToList(loads, node.Value, bcDofs, bcValues);
			                }
			                break;
		                }
	                
		                case BoundaryConditionCase.LeftDistributedFlux or BoundaryConditionCase.RightDistributedFlux or 
			                BoundaryConditionCase.FrontDistributedFlux or BoundaryConditionCase.BackDistributedFlux or 
			                BoundaryConditionCase.TopDistributedFlux or BoundaryConditionCase.BottomDistributedFlux:
		                {
			                throw new NotImplementedException("Wank!");
		                }
	                
	                }
                }
            }
        }
        var modelNeumannConditions = new StructuralBoundaryConditionSet(new List<INodalDisplacementBoundaryCondition>(), loads);

		
        model.BoundaryConditions.Add(modelNeumannConditions);
    }
    
	private static void AddNodalForceToList(List<INodalLoadBoundaryCondition> loads,
	    INode node, StructuralDof[] dofs, double[] values)
    {
	    for (int i = 0; i < dofs.Length; i++)
	    {
		    loads.Add(new NodalLoad(node, dofs[i], values[i]));
	    }
    }

    private static void AddBoundedDisplacementDofToList(List<INodalDisplacementBoundaryCondition> constraints,
		INode node, StructuralDof[] dofs, double[] values)
	{
		for (int i = 0; i < dofs.Length; i++)
		{
			constraints.Add(new NodalDisplacement(node, dofs[i], values[i]));
		}
	}
    
	
	public static void AssignConvectionDiffusionDirichletBCsToModel(Model model,
	List<(BoundaryConditionCase, ConvectionDiffusionDof[], double[][], double[])> bcs,
	double tolerance)
	{
		var constraints = new List<INodalConvectionDiffusionDirichletBoundaryCondition>();
		
		foreach (var bcSet in bcs)
		{
			var bcFace = bcSet.Item1;
			var bcDofs = bcSet.Item2;
			var faceCoords = bcSet.Item3;
			var bcValues = bcSet.Item4;
			
			foreach (var node in model.NodesDictionary)
			{
				var nodalCoords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };
				
				switch (bcFace)
				{
					case BoundaryConditionCase.LeftDirichlet or BoundaryConditionCase.RightDirichlet:
					{
						if (Math.Abs(nodalCoords[0] - faceCoords[0][0]) < tolerance)
						{
							AddBoundedConvectionDiffusionDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
					case BoundaryConditionCase.FrontDirichlet or BoundaryConditionCase.BackDirichlet:
					{
						if (Math.Abs(nodalCoords[1] - faceCoords[0][1]) < tolerance)
						{
							AddBoundedConvectionDiffusionDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
					case BoundaryConditionCase.TopDirichlet or BoundaryConditionCase.BottomDirichlet:
					{
						if (Math.Abs(nodalCoords[2] - faceCoords[0][2]) < tolerance)
						{
							AddBoundedConvectionDiffusionDofToList(constraints, node.Value, bcDofs, bcValues);
						}
						break;
					}
				}
			}
		}
		var modelDirichletConditions = new ConvectionDiffusionBoundaryConditionSet(constraints, new List<INodalConvectionDiffusionNeumannBoundaryCondition>());
		model.BoundaryConditions.Add(modelDirichletConditions);
	}
	
	private static void AddBoundedConvectionDiffusionDofToList(List<INodalConvectionDiffusionDirichletBoundaryCondition> constraints,
		INode node, ConvectionDiffusionDof[] dofs, double[] values)
	{
		constraints.Add(new NodalUnknownVariable(node, dofs[0], values[0]));
	}
	
	public static void AssignConvectionDiffusionICToModel(Model model, double initialValue)
	{
		//var modelBoundaryConditions = (model.BoundaryConditions[0] as ConvectionDiffusionBoundaryConditionSet)
		//	.EnumerateNodalBoundaryConditions().ToList() as List<INodalBoundaryCondition<IConvectionDiffusionDofType>>;
		//var dirichletNodes = (modelBoundaryConditions
		//		.Where(x => x is INodalConvectionDiffusionDirichletBoundaryCondition).ToList())
		//	.Select(x => x.Node).Distinct().ToList();
		//var freeNodes = model.NodesDictionary.Values.Except(dirichletNodes).ToList();
		//var initialConditions = new List<INodalConvectionDiffusionInitialCondition>();
		//foreach (var freeNode in freeNodes)
		//{
		//	initialConditions.Add(new NodalInitialUnknownVariable(freeNode, ConvectionDiffusionDof.UnknownVariable, initialValue));
		//}

		//model.InitialConditions.Add(new ConvectionDiffusionInitialConditionSet(initialConditions, new DomainInitialUnknownVariable[]{ })); //TODOnewnug
	}

}
