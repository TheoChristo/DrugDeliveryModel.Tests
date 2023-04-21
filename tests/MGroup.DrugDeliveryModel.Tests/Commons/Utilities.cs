using System.Collections.Generic;
using MGroup.MSolve.DataStructures;
using MGroup.NumericalAnalyzers.Logging;
using MGroup.MSolve.Discretization.Entities;
using Xunit;
using System;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.Constitutive.Structural.BoundaryConditions;
using MGroup.Constitutive.Structural.Continuum;
using MGroup.Constitutive.Structural.Transient;
using MGroup.Constitutive.Structural;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Numerics.Interpolation;
using MGroup.MSolve.Numerics.Integration.Quadratures;
using System.Xml.Linq;

namespace MGroup.DrugDeliveryModel.Tests.Commons
{


	public enum MSolveCases
	{
		p_u_OneWay_FreeTopFace,
		p_u_OneWay_FreeTopRightFaces,
		u_p_OneWay_FreeTopFace,
		u_p_OneWay_FreeTopRightFaces,
		u_P
	}

	public static class Utilities
	{


		public static Model GetParallelepipedMesh()
		{
			var nnx = 3;
			var nny = 3;
			var nnz = 3;
			var lx = 1;
			var ly = 1;
			var lz = 1;
			var maxX = 1;
			var maxY = 1;
			var maxZ = 1;
			var mesh = new OrthogonalParallelepipedHexa8MeshGenerator(nnx, nny, nnz, lx, ly, lz, maxX, maxY, maxZ);
			var nodes = mesh.Nodes;

			//var mesh = new ComsolMeshReader(fileName);
			//var nodes = mesh.NodesDictionary;

			var model = new Model();
			model.SubdomainsDictionary[0] = new Subdomain(id: 0);

			foreach (var node in nodes.Values)
			{
				model.NodesDictionary.Add(node.ID, node);
			}

			var elasticMaterial1 = new ElasticMaterial3D(youngModulus: 1, poissonRatio: 0.3);
			var elementFactory =
				new ContinuumElement3DFactory(elasticMaterial1, new TransientAnalysisProperties(1, 0, 0));

			//var domains = new Dictionary<int, double[]>(2);
			foreach (var elementConnectivity in mesh.ElementConnectivity)
			{
				var domainId = elementConnectivity.Value.Item3;
				var element =
					elementFactory.CreateElement(elementConnectivity.Value.Item1, elementConnectivity.Value.Item2);
				model.ElementsDictionary.Add(elementConnectivity.Key, element);
				model.SubdomainsDictionary[0].Elements.Add(element);
			}

			var faceXYNodes = new List<INode>();
			var faceXZNodes = new List<INode>();
			var faceYZNodes = new List<INode>();

			foreach (var node in model.NodesDictionary.Values)
			{
				if (Math.Abs(0 - node.Z) < 1E-9) faceXYNodes.Add(node);
				if (Math.Abs(0 - node.Y) < 1E-9) faceXZNodes.Add(node);
				if (Math.Abs(0 - node.X) < 1E-9) faceYZNodes.Add(node);
			}

			var constraints = new List<INodalDisplacementBoundaryCondition>();
			foreach (var node in faceXYNodes)
			{
				constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationZ, amount: 0d));
			}

			foreach (var node in faceXZNodes)
			{
				constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationY, amount: 0d));
			}

			foreach (var node in faceYZNodes)
			{
				constraints.Add(new NodalDisplacement(node, StructuralDof.TranslationX, amount: 0d));
			}

			INode maxDistanceNode = null;
			double currentMaxDistance = 0;
			foreach (INode node in model.NodesDictionary.Values)
			{
				double distance = Math.Sqrt(Math.Pow(node.X, 2) + Math.Pow(node.Y, 2) + Math.Pow(node.Z, 2));
				if (distance > currentMaxDistance)
				{
					currentMaxDistance = distance;
					maxDistanceNode = node;
				}
			}

			var loads = new List<INodalLoadBoundaryCondition>();

			loads.Add(new NodalLoad
			(
				maxDistanceNode,
				StructuralDof.TranslationZ,
				amount: 0.00000001
			));

			model.BoundaryConditions.Add(new StructuralBoundaryConditionSet(constraints, loads));

			return model;
		}

		public static bool AreDisplacementsSame(IReadOnlyList<double[]> expectedDisplacements,
			TotalDisplacementsPerIterationLog computedDisplacements, double tolerance)
		{
			var comparer = new ValueComparer(tolerance);
			for (var iter = 0; iter < expectedDisplacements.Count; ++iter)
			{
				for (var i = 0; i < expectedDisplacements[iter].Length; ++i)
				{
					var expected = expectedDisplacements[iter][i];
					(var node, var dof) = computedDisplacements.WatchDofs[i];
					var computed = computedDisplacements.GetTotalDisplacement(iter, node, dof);

					if (!comparer.AreEqual(expected, computed)) return false;
				}
			}

			return true;
		}

		public static bool AreDisplacementsSame(IReadOnlyList<double[]> expectedDisplacements,
			IncrementalDisplacementsLog computedDisplacements, double tolerance)
		{
			var comparer = new ValueComparer(tolerance);
			for (var iter = 0; iter < expectedDisplacements.Count; ++iter)
			{
				for (var i = 0; i < expectedDisplacements[iter].Length; ++i)
				{
					var expected = expectedDisplacements[iter][i];
					(var node, var dof) = computedDisplacements.WatchDofs[i];
					var computed = computedDisplacements.GetTotalDisplacement(iter, node, dof);

					if (!comparer.AreEqual(expected, computed)) return false;
				}
			}

			return true;
		}

		public static bool AreDisplacementsSame(double[] expectedDisplacements, double[] computedDisplacements,
			double tolerance)
		{
			var comparer = new ValueComparer(tolerance);

			for (var i = 0; i < expectedDisplacements.Length; i++)
			{
				if (!comparer.AreEqual(expectedDisplacements[i], computedDisplacements[i])) return false;
			}

			return true;
		}

		public static void CheckModelSubdomains(Dictionary<int, int[]> expectedSubdomains, Model model)
		{
			for (var i = 0; i < expectedSubdomains.Count; i++)
			{
				var subdomainElements = model.SubdomainsDictionary[i].Elements;
				Assert.Equal(expectedSubdomains[i].Length, model.SubdomainsDictionary[i].Elements.Count);
				for (var j = 0; j < expectedSubdomains[i].Length; j++)
				{
					Assert.Equal(expectedSubdomains[i][j], subdomainElements[j].ID);
				}
			}
		}

		public static bool AreTensorsEqual(IReadOnlyList<double[]> tensors1, IReadOnlyList<double[]> tensors2,
			double tolerance)
		{
			if (tensors1.Count != tensors2.Count) return false;
			for (int i = 0; i < tensors1.Count; ++i)
			{
				if (tensors1[i].Length != tensors2[i].Length) return false;
				for (int j = 0; j < tensors1[i].Length; ++j)
				{
					if (!AreValuesEqual(tensors1[i][j], tensors2[i][j], tolerance)) return false;
				}
			}

			return true;
		}

		public static bool AreValuesEqual(double value1, double value2, double tolerance)
		{
			if (Math.Abs(value2) <= tolerance) // Can't divide with expected ~= 0. 
			{
				if (Math.Abs(value1) <= tolerance) return true;
				else return false;
			}
			else return (Math.Abs(1.0 - value1 / value2) < tolerance) ? true : false;
		}

		public static int FindNodeIdFromNodalCoordinates(Dictionary<int, Node> nodes, double[] nodalCoords,
			double tolerrance)
		{
			var id = -1;
			var currentDistance = 1000d;
			foreach (var node in nodes.Values)
			{
				double distance = Math.Sqrt(Math.Pow(node.X - nodalCoords[0], 2) +
				                            Math.Pow(node.Y - nodalCoords[1], 2) +
				                            Math.Pow(node.Z - nodalCoords[2], 2));
				if (distance < currentDistance)
				{
					currentDistance = distance;
					id = node.ID;
				}
			}

			Console.WriteLine("Node found with ID: " + id + " and Coordinates: " + nodes[id].X + " " + nodes[id].Y +
			                  " " +
			                  nodes[id].Z);
			return id;
		}

		public static int FindElementIdFromGaussPointCoordinates(Model model, double[] gpCoords, double tolerrance)
		{
			var id = -1;
			var cuurentDistance = 1000d;
			foreach (var element in model.ElementsDictionary)
			{
				var elementGpCoords = GetGaussPointsCoordinates(0, element.Value);
                //= ((ConvectionDiffusionElement3D)element.Value).GetGaussPointsCoordinates(0);
                double distance = Math.Sqrt(Math.Pow(elementGpCoords[0] - gpCoords[0], 2) +
				                            Math.Pow(elementGpCoords[1] - gpCoords[1], 2) +
				                            Math.Pow(elementGpCoords[2] - gpCoords[2], 2));
				if (distance < cuurentDistance)
				{
					cuurentDistance = distance;
					id = element.Key;
				}
			}

			var foundGpCoords = GetGaussPointsCoordinates(0, model.ElementsDictionary[id]);
               // ((ConvectionDiffusionElement3D)model.ElementsDictionary[id]).GetGaussPointsCoordinates(0);
			Console.WriteLine("GP with coordinates: " + foundGpCoords[0] + " " + foundGpCoords[1] + " " +
			                  foundGpCoords[2] +
			                  " is in element with id: " + id);
			return id;
		}

		public static double[] GetGaussPointsCoordinates(int gpNo, IElementType element)
		{
			var interpolation = InterpolationTet4.UniqueInstance;
			var quadrature = TetrahedronQuadrature.Order1Point1;
			var Nodes = element.Nodes;
			var shapeFunctionValues = interpolation.EvaluateFunctionsAt(quadrature.IntegrationPoints[gpNo]);
			double[] gpCoordinates = new double[3]; //{ dphi_dksi, dphi_dheta, dphi_dzeta}
			for (int i1 = 0; i1 < shapeFunctionValues.Length; i1++)
			{
				gpCoordinates[0] += shapeFunctionValues[i1] * Nodes[i1].X;
				gpCoordinates[1] += shapeFunctionValues[i1] * Nodes[i1].Y;
				gpCoordinates[2] += shapeFunctionValues[i1] * Nodes[i1].Z;
			}

			return gpCoordinates;
		}

		public static List<INode> FindElementNodesFromElementId(Model model, int elementId)
		{
			var element = model.ElementsDictionary[elementId];
			var nodes = new List<INode>();
			for (int i = 0; i < element.Nodes.Count; i++)
			{
				nodes.Add(element.Nodes[i]);
				Console.WriteLine("Node ID: " + nodes[i].ID + " Node X: " + nodes[i].X + " Node Y: " + nodes[i].Y +
				                  " Node Z: " + nodes[i].Z);
			}

			return nodes;
		}

		public static int FindRandomInternalNode(Dictionary<int, Node> nodes, double minX, double maxX,
			double minY, double maxY,
			double minZ, double maxZ)
		{
			var id = -1;

			foreach (var node in nodes)
			{
				var coords = new double[] { node.Value.X, node.Value.Y, node.Value.Z };
				if (coords[0] > minX && coords[0] < maxX &&
				    coords[1] > minY && coords[1] < maxY &&
				    coords[2] > minZ && coords[2] < maxZ)
				{
					id = node.Key;
					break;
				}

				id = node.Key;
			}

			Console.WriteLine("Internal Node found with ID: " + id + " and Coordinates: " + nodes[id].X + " " +
			                  nodes[id].Y + " " +
			                  nodes[id].Z);
			return id;
		}
	}
}