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
using MGroup.Constitutive.Structural;
using MGroup.DrugDeliveryModel.Tests.Commons;
using Xunit;
using MGroup.Constitutive.ConvectionDiffusion;
using MGroup.FEM.ConvectionDiffusion.Isoparametric;
using MGroup.FEM.Structural.Continuum;
using TriangleNet.Meshing.Algorithm;

namespace MGroup.DrugDeliveryModel.Tests.Integration
{
    public class Coupled7and9eqsSolutionex6_1
    {
        const double Sc = 0.1;

        private const double timeStep = 1.0e-5; // in sec
        const double totalTime = 4.0e-4; // in sec
        static int incrementsPertimeStep = 1;


        // strucutral model Loads
        static StructuralDof loadedDof = StructuralDof.TranslationZ;
        static double load_value = 0; //[kN]


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


        static double pressureMonitorNodeX = 0.05;
        static double pressureMonitorNodeY = 0.05;
        static double pressureMonitorNodeZ = 0.05;

        static double[] pressureMonitorNode = new double[]
            { pressureMonitorNodeX, pressureMonitorNodeY, pressureMonitorNodeZ };

        private static int pressureMonitorID;


        //Mid node top face (z = zmax)
        static double divPMonitorGPX = 0.05;
        static double divPMonitorGPY = 0.05;
        static double divPMonitorGPZ = 0.1;

        //Mid node right face (x = xmax)
        //static double divPMonitorGPX = 0.1;
        //static double divPMonitorGPY = 0.05;
        //static double divPMonitorGPZ = 0.05;

        static double[] divPMonitorGP = new double[] { divPMonitorGPX, divPMonitorGPY, divPMonitorGPZ };

        private static int divPMonitorID;

        static StructuralDof eq9dofTypeToMonitor = StructuralDof.TranslationX;

        //Darcy model properties

        //static double k_th = 7.52e-14; // m2/(KPa sec)
        static double k_th = 7.52e-10; // m2/(KPa sec)
        //static double k_th = 7.5231E-8 ; // m2/(KPa sec)
        static double Sv = 7e+3; // 1/(m)
        static double Lp = 2.7e-12; // m/(KPa sec)
        static double LplSvl = 0; // 1/(KPa sec)
        static double pv = 4; // kPa
        static double pl = 0d; // KPa
        static double div_vs = 1e-6; // 1/(sec)

        //static int nodeIdToMonitor = 36;

        //The above coordinates are in the reference configuration.
        //The actual coordinates will be computed in the analysis.

        static double structuralMonitorNodeX = 0.05;
        static double structuralMonitorNodeY = 0.05;
        static double structuralMonitorNodeZ = 0.1;

        static double[] structuralMonitorNode = new double[]
            { structuralMonitorNodeX, structuralMonitorNodeY, structuralMonitorNodeZ };

        private static int structuralMonitorID;
        static ConvectionDiffusionDof eq7n8dofTypeToMonitor = ConvectionDiffusionDof.UnknownVariable;

        //private static double boundaryValueAllBoundaries = pv;
        //private static double initialCondition = pv;
        private static double boundaryValueAllBoundaries = 0;
        private static double initialCondition = 0;
        //private static double boundaryValueAllBoundaries = 1;
        //private static double initialCondition = 1;

        static double modelMinX = 0;
        static double modelMaxX = 0.1;
        static double modelMinY = 0;
        static double modelMaxY = 0.1;
        static double modelMinZ = 0;
        static double modelMaxZ = 0.1;

        static double[] expextedPnode = new double[] { 0.0548916067892445, 0.0548916069686221, 0.0548916072735233, 0.0548916077128598, 0.054891608291035, 0.0548916090101852, 0.0548916098713005, 0.054891610874785, 0.0548916120207372, 0.05489161330909, 0.0548916147396809, 0.0548916163122873, 0.054891618026644, 0.0548916198824525, 0.0548916218793851, 0.0548916240170879, 0.0548916262951818, 0.0548916287132634, 0.0548916312709056, 0.054891633967658, 0.0548916368030475, 0.0548916397765784, 0.0548916428877328, 0.0548916461359713, 0.054891649520733, 0.0548916530414361, 0.0548916566974781, 0.0548916604882365, 0.0548916644130689, 0.0548916684713138, 0.0548916726622905, 0.0548916769853001, 0.0548916814396256, 0.0548916860245323, 0.0548916907392685, 0.0548916955830657, 0.0548917005551392, 0.0548917056546886, 0.0548917108808979, 0.0548917162329363, 0.0548917217099588, 0.0548917273111059, 0.0548917330355049, 0.0548917388822698, 0.0548917448505019, 0.0548917509392903, 0.0548917571477121, 0.0548917634748329, 0.0548917699197076, 0.05489177648138, 0.054891783158884, 0.0548917899512434, 0.0548917968574727, 0.054891803876577, 0.0548918110075528, 0.0548918182493882, 0.0548918256010632, 0.0548918330615499, 0.054891840629813, 0.0548918483048102, 0.054891856085492, 0.0548918639708028, 0.0548918719596804, 0.0548918800510564, 0.0548918882438571, 0.0548918965370028, 0.0548919049294088, 0.0548919134199849, 0.0548919220076365, 0.054891930691264, 0.0548919394697633, 0.054891948342026, 0.0548919573069397, 0.0548919663633878, 0.0548919755102499, 0.0548919847464018, 0.054891994070716, 0.0548920034820613, 0.0548920129793033, 0.0548920225613042, 0.0548920322269233, 0.0548920419750168, 0.0548920518044379, 0.0548920617140371, 0.054892071702662, 0.0548920817691577, 0.0548920919123667, 0.0548921021311286, 0.0548921124242812, 0.0548921227906593, 0.0548921332290956, 0.0548921437384208, 0.0548921543174631, 0.0548921649650486, 0.0548921756800013, 0.0548921864611435, 0.0548921973072953, 0.0548922082172748, 0.0548922191898987 };
        static double[] expectedDisplacementXNode = new double[] { 2.583948896393967e-12, 9.038478059870692e-12, 1.9997828876538827e-11, 3.576606186015938e-11, 5.6478045381061846e-11, 8.217990640777271e-11, 1.1286925116404781e-10, 1.48515292994697e-10, 1.8906894775777576e-10, 2.344679227153957e-10, 2.8463930980201074e-10, 3.3950093763238794e-10, 3.9896210887231685e-10, 4.62924035970249e-10, 5.312801315709562e-10, 6.039162316481332e-10, 6.80710790275793e-10, 7.615350654186316e-10, 8.46253305344712e-10, 9.347229403852783e-10, 1.026794782310252e-09, 1.1223132323689511e-09, 1.2211164984227483e-09, 1.3230368212796684e-09, 1.427900710182029e-09, 1.5355291873063727e-09, 1.6457380410906255e-09, 1.7583380881621775e-09, 1.8731354436323717e-09, 1.9899317994963443e-09, 2.108524710864787e-09, 2.228707889739846e-09, 2.350271506038031e-09, 2.473002495542235e-09, 2.596684874458469e-09, 2.721100060237029e-09, 2.846027198305762e-09, 2.9712434943495156e-09, 3.0965245517569947e-09, 3.2216447138439556e-09, 3.3463774104477413e-09, 3.470495508478038e-09, 3.5937716659948286e-09, 3.7159786893694923e-09, 3.836889893080333e-09, 3.9562794616741106e-09, 4.073922813421889e-09, 4.1895969651778e-09, 4.303080897948238e-09, 4.414155922660666e-09, 4.522606045620081e-09, 4.628218333119754e-09, 4.7307832746828055e-09, 4.830095144384584e-09, 4.925952359711221e-09, 5.018157837398237e-09, 5.1065193456886326e-09, 5.190849852444281e-09, 5.2709678685424555e-09, 5.346697785982421e-09, 5.417870210127226e-09, 5.484322285504869e-09, 5.545898014590985e-09, 5.602448568995417e-09, 5.6538325924760356e-09, 5.699916495208093e-09, 5.740574738739476e-09, 5.775690111064371e-09, 5.805153991258223e-09, 5.828866603114857e-09, 5.846737257248619e-09, 5.8586845811156305e-09, 5.86463673643301e-09, 5.864531623479105e-09, 5.8583170717686785e-09, 5.8459510166120755e-09, 5.827401661084154e-09, 5.802647622938445e-09, 5.771678066027287e-09, 5.734492815793515e-09, 5.6911024584314725e-09, 5.641528423324461e-09, 5.585803048393057e-09, 5.523969628005016e-09, 5.45608244312539e-09, 5.38220677340493e-09, 5.302418890931225e-09, 5.216806035393306e-09, 5.125466370432959e-09, 5.0285089209864925e-09, 4.9260534914499645e-09, 4.818230564522131e-09, 4.705181180618075e-09, 4.587056797765477e-09, 4.464019131938743e-09, 4.336239977804984e-09, 4.203901009894459e-09, 4.067193564237599e-09, 3.926318400542259e-09 };
        static double[] expectedDisplacementYNode = new double[] { -1.3698329043362804e-12, -4.795188920330931e-12, -1.0620222350026239e-11, -1.9018885043680948e-11, -3.008060024757297e-11, -4.385309330471974e-11, -6.036379225224811e-11, -7.9630507279698e-11, -1.0166674478850303e-10, -1.264843346141752e-10, -1.5409471008517582e-10, -1.8450951102192918e-10, -2.1774084493375336e-10, -2.5380137421786285e-10, -2.927043134325794e-10, -3.3446337885855286e-10, -3.7909271162149724e-10, -4.2660678520862817e-10, -4.770203029798862e-10, -5.303480886765171e-10, -5.866049716324947e-10, -6.458056677564068e-10, -7.079646570372229e-10, -7.730960581718182e-10, -8.412135008429377e-10, -9.123299961450263e-10, -9.864578056438292e-10, -1.0636083095436754e-09, -1.1437918744450795e-09, -1.227017721171506e-09, -1.3132937931402672e-09, -1.4026266257570366e-09, -1.4950212173102995e-09, -1.590480901832853e-09, -1.6890072243991743e-09, -1.7905998193142766e-09, -1.895256291646697e-09, -2.002972102541809e-09, -2.113740458746077e-09, -2.2275522067567687e-09, -2.344395731995211e-09, -2.464256863389009e-09, -2.587118783728242e-09, -2.7129619461428902e-09, -2.841763997031417e-09, -2.9734997057427804e-09, -3.108140901297363e-09, -3.2456564164061983e-09, -3.386012039023381e-09, -3.529170471639303e-09, -3.6750912985033163e-09, -3.823730960922399e-09, -3.975042740773482e-09, -4.1289767523224206e-09, -4.285479942423409e-09, -4.444496099135725e-09, -4.605965868771698e-09, -4.769826781360345e-09, -4.936013284470746e-09, -5.104456785319248e-09, -5.275085701052895e-09, -5.447825517066802e-09, -5.622598853189558e-09, -5.799325537537754e-09, -5.977922687813242e-09, -6.158304799785664e-09, -6.340383842681864e-09, -6.524069361172879e-09, -6.709268583620803e-09, -6.895886536228076e-09, -7.0838261627074264e-09, -7.272988449062747e-09, -7.463272553057346e-09, -7.654575937922837e-09, -7.846794509838532e-09, -8.039822758704342e-09, -8.233553901704748e-09, -8.427880029152443e-09, -8.622692252090485e-09, -8.817880851113111e-09, -9.0133354258664e-09, -9.208945044675212e-09, -9.404598393745732e-09, -9.60018392538012e-09, -9.79559000464725e-09, -9.990705053953182e-09, -1.0185417694952738e-08, -1.0379616887255664e-08, -1.05731920633835e-08, -1.0766033259440904e-08, -1.0958031240982994e-08, -1.1149077623565406e-08, -1.133906498748512e-08, -1.1527886986233661e-08, -1.17154384482025e-08, -1.1901615471205765e-08, -1.2086315509402777e-08, -1.226943745223102e-08, -1.2450881694983224e-08 };
        static double[] expectedDisplacementZNode = new double[] { 1.1538114693538292e-10, 4.0382617091610815e-10, 8.941630891320444e-10, 1.600786645864557e-09, 2.530868808379729e-09, 3.687963880343075e-09, 5.073810964472423e-09, 6.689235115793462e-09, 8.534597844006441e-09, 1.0610022288214e-08, 1.2915505722745617e-08, 1.5450975722108147e-08, 1.8216318148538852e-08, 2.121139104376313e-08, 2.4436031465835773e-08, 2.7890058791704853e-08, 3.1573276246080193e-08, 3.5485471537200026e-08, 3.96264170401366e-08, 4.3995869748313746e-08, 4.859357110395887e-08, 5.341924676324317e-08, 5.84726063244012e-08, 6.37533430334046e-08, 6.926113347491235e-08, 7.499563725281604e-08, 8.095649666299926e-08, 8.714333636009884e-08, 9.355576301963832e-08, 1.0019336499671273e-07, 1.0705571198231357e-07, 1.1414235465833949e-07, 1.2145282435232057e-07, 1.289866326928858e-07, 1.367432712670003e-07, 1.447222112799972e-07, 1.5292290321943854e-07, 1.6134477652382606e-07, 1.6998723925719032e-07, 1.7884967779057445e-07, 1.879314564914261e-07, 1.972319174218947e-07, 2.067503800470295e-07, 2.1648614095384545e-07, 2.2643847358222e-07, 2.366066279685509e-07, 2.469898305031e-07, 2.575872837019027e-07, 2.683981659941166e-07, 2.7942163152564004e-07, 2.906568099798067e-07, 3.021028064159264e-07, 3.137587011264084e-07, 3.256235495131632e-07, 3.37696381983944e-07, 3.49976203869244e-07, 3.624619953603232e-07, 3.7515271146889257e-07, 3.8804728200893913e-07, 4.01144611601122e-07, 4.144435797001271e-07, 4.279430406453098e-07, 4.416418237349054e-07, 4.555387333240309e-07, 4.696325489466517e-07, 4.839220254616258e-07, 4.98405893222881e-07, 5.130828582737285e-07, 5.279516025652507e-07, 5.430107841986459e-07, 5.582590376913557e-07, 5.736949742667344e-07, 5.893171821669678e-07, 6.051242269888846e-07, 6.21114652042247e-07, 6.372869787300417e-07, 6.536397069502439e-07, 6.701713155184611e-07, 6.868802626108086e-07, 7.037649862263116e-07, 7.208239046680727e-07, 7.380554170423941e-07, 7.554579037749813e-07, 7.730297271433125e-07, 7.907692318242055e-07, 8.086747454555604e-07, 8.267445792112177e-07, 8.449770283878208e-07, 8.633703730025348e-07, 8.819228784004256e-07, 9.006327958702756e-07, 9.194983632675656e-07, 9.385178056433356e-07, 9.57689335877584e-07, 9.770111553158629e-07, 9.964814544076826e-07, 1.0160984133453336e-06, 1.0358602027016948e-06, 1.0557649840656092e-06 };

        public Coupled7and9eqsSolutionex6_1()
        {
            IsoparametricJacobian3D.DeterminantTolerance = 1e-20;
        }

        [Theory]
        //[InlineData("../../../DataFiles/workingTetMesh4886.mphtxt")]
        //[InlineData("../../../DataFiles/chipMelter2M.mphtxt")]
        //[InlineData("../../../DataFiles/MeshCyprusTM.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh4886_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh648_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingHexaMesh64_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingHexaMesh8_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh155.mphtxt")]
        //[InlineData("../../../DataFiles/workingTetMesh_1Domain.mphtxt")]
        //[InlineData("../../../DataFiles/workingHexaMesh216_1Domain.mphtxt")]
        [InlineData("../../../DataFiles/workingTetMesh2185_1Domain.mphtxt")]


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
            structuralMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, structuralMonitorNode, 1e-2);
            pressureMonitorID = Utilities.FindNodeIdFromNodalCoordinates(comsolReader.NodesDictionary, pressureMonitorNode, 1e-2);

            //pressureMonitorID = Utilities.FindRandomInternalNode(comsolReader.NodesDictionary, modelMinX, modelMaxX,
            //modelMinY, modelMaxY, modelMinZ, modelMaxZ);
            var eq78Model = new Eq78ModelProviderForStaggeredSolutionex6_1(comsolReader, k_th, Lp, Sv, pv, LplSvl, pl,
                velocityDivergenceAtElementGaussPoints,

                boundaryValueAllBoundaries, initialCondition, pressureMonitorID, eq7n8dofTypeToMonitor, modelMinX,
                modelMaxX, modelMinY, modelMaxY, modelMinZ, modelMaxZ);

            var eq9Model = new Eq9ModelProviderForStaggeredSolutionex6_1(comsolReader, Sc, miNormal, kappaNormal, miTumor,
                kappaTumor, timeStep, totalTime, lambda, pressureTensorDivergenceAtElementGaussPoints,
                structuralMonitorID, eq9dofTypeToMonitor, loadedDof, load_value, modelMinX, modelMaxX, modelMinY,
                modelMaxY, modelMinZ, modelMaxZ);

            var equationModel = new Coupled7and9eqsModelex6_1(eq78Model, eq9Model, comsolReader, lambda,
                pressureTensorDivergenceAtElementGaussPoints, velocityDivergenceAtElementGaussPoints, timeStep,
                totalTime, incrementsPertimeStep);



            #region loggin

            var p_i = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsX = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsY = new double[(int)(totalTime / timeStep)];
            double[] structuralResultsZ = new double[(int)(totalTime / timeStep)];
            var displacements = new List<double[]>();
            displacements.Add(structuralResultsX);
            displacements.Add(structuralResultsY);
            displacements.Add(structuralResultsZ);

            double[] modelMaxVelDivOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dxOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dyOverTime = new double[(int)(totalTime / timeStep)];
            //double[] modelMax_dP_dzOverTime = new double[(int)(totalTime / timeStep)];
            //var dp_dxi = new List<double[]>();
            //dp_dxi.Add(modelMax_dP_dxOverTime);
            //dp_dxi.Add(modelMax_dP_dyOverTime);
            //dp_dxi.Add(modelMax_dP_dzOverTime);
            #endregion


            var staggeredAnalyzer = new StepwiseStaggeredAnalyzer(equationModel.ParentAnalyzers,
                equationModel.ParentSolvers, equationModel.CreateModel, maxStaggeredSteps: 200, tolerance: 0.001);
            for (currentTimeStep = 0; currentTimeStep < totalTime / timeStep; currentTimeStep++)
            {
                equationModel.CurrentTimeStep = currentTimeStep;
                equationModel.CreateModelFirstTime(equationModel.ParentAnalyzers, equationModel.ParentSolvers);
                staggeredAnalyzer.SolveCurrentStep();

                #region logging

                var allValues = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues
                    .Select(x => x.val).ToArray();

                var divPMonitorId = Utilities.FindElementIdFromGaussPointCoordinates(equationModel.model[0], divPMonitorGP, 1e-1);

                p_i[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[0].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[0].GetNode(pressureMonitorID), ConvectionDiffusionDof.UnknownVariable];
                structuralResultsX[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationX];
                structuralResultsY[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationY];
                structuralResultsZ[currentTimeStep] = ((DOFSLog)equationModel.ParentAnalyzers[1].ChildAnalyzer.Logs[0]).DOFValues[equationModel.model[1].GetNode(structuralMonitorID), StructuralDof.TranslationZ];


                //modelMax_dP_dxOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .xcoeff_OverTimeAtGp1[0];
                //modelMax_dP_dyOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .ycoeff_OverTimeAtGp1[0];
                //modelMax_dP_dzOverTime[currentTimeStep] =
                //    ((ConvectionDiffusionElement3D)equationModel.model[0].ElementsDictionary[divPMonitorId])
                //    .zcoeff_OverTimeAtGp1[0];



                modelMaxVelDivOverTime[currentTimeStep] = velocityDivergenceAtElementGaussPoints
                    .Select(x => Math.Abs(x.Value[0])).ToArray().Max();

                //modelMax_dP_dxOverTime[currentTimeStep] = pressureTensorDivergenceAtElementGaussPoints
                //    .Select(x => Math.Abs(x.Value[0][0])).ToArray().Max();


                if (Solution.ContainsKey(currentTimeStep))
                {
                    Solution[currentTimeStep] = allValues;
                    Console.WriteLine($"Time step: {timeStep}");
                    Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[timeStep])}");
                }
                else
                {
                    Solution.Add(currentTimeStep, allValues);
                }

                #endregion

                (equationModel.ParentAnalyzers[0] as NewmarkDynamicAnalyzer).AdvanceStep();
                (equationModel.ParentAnalyzers[1] as NewmarkDynamicAnalyzer).AdvanceStep();

                for (int j = 0; j < equationModel.ParentAnalyzers.Length; j++)
                {
                    equationModel.AnalyzerStates[j] = equationModel.ParentAnalyzers[j].CreateState();
                    equationModel.NLAnalyzerStates[j] = equationModel.NLAnalyzers[j].CreateState();
                }

                Console.WriteLine($"Displacement vector: {string.Join(", ", Solution[currentTimeStep])}");
            }

            double pr = 100;
            double Fval = 100;
            double Fval_e = 0;
            //eq9model.load_value;


/*            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsX,
                $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{1}_.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsY,
                $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{2}_.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsZ,
                $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_u_{3}_.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(modelMaxVelDivOverTime,
                $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_modelMaxVelDivOverTime.txt");
            (new MGroup.LinearAlgebra.Output.Array1DWriter()).WriteToFile(structuralResultsZ,
                $@"C:\Users\acivi\Documents\atuxaia\develop yperelastic withh BIO TEAM\VALIDATION EQs1\staggered_ex5\pressure_{pr}_F_{Fval}_e{Fval_e}_LOGGEDval_modelMax_dP_dx{3}OverTime.txt");
*/


            var path = "../../../StaggeredSolutionPresDynamex6/6.1/dp_dxi_mslv.csv";
            //CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(dp_dxi), "../../../StaggeredSolutionPresDynamex6/6.1/dp_dxi_GP_mslv.csv");
            CSVExporter.ExportMatrixToCSV(CSVExporter.ConverVectorsTo2DArray(displacements), "../../../StaggeredSolutionPresDynamex6/6.1/displacements_nodes_mslv.csv");
            CSVExporter.ExportVectorToCSV(p_i, "../../../StaggeredSolutionPresDynamex6/pi_nodes_mslv.csv");

            Assert.True(CompareResults(p_i, displacements));

        }

        private bool CompareResults(double[] pressure, List<double[]> displacements)
        {
            bool ret = true;
            double tolerance = 1e-7;

            for (var i = 0; i < expextedPnode.Length; i++)
            {
                double pError = Math.Abs(pressure[i] - expextedPnode[i]);
                double dispXError = Math.Abs(displacements[0][i] - expectedDisplacementXNode[i]);
                double dispYError = Math.Abs(displacements[1][i] - expectedDisplacementYNode[i]);
                double dispZError = Math.Abs(displacements[2][i] - expectedDisplacementZNode[i]);

                if (pError > tolerance || dispXError > tolerance || dispYError > tolerance || dispZError > tolerance)
                {
                    Console.WriteLine("Wrong result on step " + i);
                    ret = false;
                    break;
                }
            }

            return ret;
        }




    }
}
