package swmutsel.model;

import swmutsel.model.parameters.Fitness;

/**
 * 
 * https://compbio.soe.ucsc.edu/dirichlets/recode3.20comp
 * 
 * NB original alphas given in alphabetical order by character code,
 * but here using alphabetical order by full length name
 *
 * @author Christopher Monit c.monit.12@ucl.ac.uk
 */
public class SAMRecode320CompPenalty extends DirichletMixPenalty {
    private static final double[][] ALPHAS = new double[][]{
        //A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
        {0.668781,0.636552,0.979541,1.08101,0.072626,0.555784,0.930936,0.982387,0.366701,0.0586812,0.185073,0.849019,0.0648062,0.1235,0.519704,1.04987,0.60931,0.0676805,0.201772,0.129148},
        {0.880009,0.581516,0.368825,0.303744,0.197393,0.439775,0.497357,0.391988,0.348159,0.936094,1.42844,0.527206,0.449302,0.657492,0.384576,0.622624,0.747064,0.235012,0.606928,1.08007},
        {0.153384,0.0252186,0.0190716,0.0073824,0.0520756,0.0274578,0.0158439,0.025533,0.0185789,0.845361,2.42256,0.0282996,0.424296,0.428964,0.0313429,0.028514,0.0519217,0.0279653,0.0664755,0.522946},
        {0.794007,0.933444,0.308302,0.624236,0.0130659,0.953727,1.85769,0.115707,0.123504,0.22099,0.341371,1.52605,0.111114,0.0290214,0.263545,0.554741,0.604551,0.00823516,0.0420054,0.396451},
        {0.740015,0.0301022,0.0300644,0.0213261,0.187165,0.0240143,0.0456854,0.0633687,0.0170331,1.06684,0.733524,0.0380614,0.138456,0.118944,0.0718692,0.0862989,0.367283,0.0113856,0.045079,1.70735},
        {0.15978,1.44604,0.171028,0.0505181,0.0261585,0.290476,0.125524,0.102549,0.157461,0.0795041,0.189383,1.26261,0.0550608,0.0350331,0.0844169,0.129158,0.138972,0.0159134,0.0637679,0.0851144},
        {0.308434,0.104878,0.213663,1.69731,0.0137217,0.320245,1.92422,0.162357,0.07232,0.0487895,0.0809074,0.236135,0.0286236,0.0361113,0.181631,0.218398,0.141668,0.0141705,0.0453433,0.0747719},
        {0.00260287,0.00459583,0.00338172,0.00631292,1.00E-05,0.00173574,0.00445502,1.04E-05,1.03E-05,0.000913052,0.00353485,0.0029241,0.00105128,0.00274753,1.04E-05,0.00274255,0.00247625,1.02E-05,0.00288489,0.00175366},
        {1.61043,0.0330646,0.0357917,0.0378292,0.15522,0.0367156,0.0498243,0.529136,0.0217524,0.040597,0.100193,0.0413396,0.0509779,0.0406484,0.0931204,0.529587,0.196607,0.00909518,0.0329275,0.230878},
        {0.15525,0.0477247,0.127437,0.0857138,0.0136827,0.0345064,0.0508316,3.10555,0.027169,0.0140491,0.0257501,0.0654038,0.00901049,0.0151451,0.0423873,0.12452,0.0341196,0.00930115,0.0187464,0.0230637},
        {0.225739,0.0751969,0.230216,0.101072,0.0684326,0.0694702,0.0813791,0.0915218,0.0336807,0.0833114,0.0731542,0.0931673,0.0419314,0.0298832,0.087446,1.13857,1.63158,0.00912576,0.0311963,0.179083},
        {1.34E-06,0.0158968,0.00678027,0.00743068,0.0166656,0.0100244,1.35E-06,0.00407061,0.00714122,1.98221,0.816669,0.017522,0.114773,0.169076,0.0106392,0.00879658,0.0399043,0.0150671,0.0517801,1.81043},
        {0.063525,0.0312317,0.752755,1.09265,0.0288391,0.0445321,0.0959581,0.216914,0.0730986,0.0207325,0.0315107,0.08719,0.00960293,0.00965196,0.059914,0.287327,0.116896,0.00701663,0.0305859,0.0249446},
        {0.294281,0.101856,0.0601581,0.12293,0.019271,0.0809085,0.162747,0.145029,0.0412349,0.0815261,0.151631,0.157594,0.021412,0.0373667,3.6966,0.23533,0.135424,0.00900473,0.0321389,0.140532},
        {0.0844832,0.0598024,0.0794234,0.0411628,0.0584945,0.0305672,0.045719,0.0590839,0.250253,0.0675757,0.168617,0.0562614,0.0439737,0.847822,0.028301,0.0798202,0.0585385,0.227395,1.30336,0.0858243},
        {0.0634034,0.00802297,1.36E-06,1.34E-06,0.0246167,0.00704479,0.00389272,0.00796028,1.35E-06,0.233395,0.541933,1.34E-06,0.101309,1.12953,0.027467,0.0248977,0.0276933,0.183309,0.516892,0.185467},
        {0.123696,0.219233,0.098734,0.0386434,0.0454619,0.554971,0.351847,0.0439442,0.223229,0.01302,0.19001,0.148699,0.120964,0.0560181,5.90E-06,0.0453885,0.0564686,0.0410248,0.0800036,0.0614792},
        {0.0212037,0.00677712,0.0105326,0.00745627,3.18769,0.00505819,0.00382411,0.0126233,1.34E-06,0.00724293,0.00785563,0.00522979,0.00489521,0.00691924,0.0136265,0.0251744,0.0235516,0.00187667,0.0038893,0.0371462},
        {0.0229376,0.0342359,0.00899225,0.00959934,0.00427768,0.00600945,0.013608,0.0227654,0.0157344,0.0226783,0.0803491,0.011561,0.0154283,0.182277,0.00980608,0.0216842,0.0189306,1.83914,0.154565,0.0223176},
        {2.17E-06,0.0252951,0.0839371,0.0198496,2.16E-06,0.0438887,2.18E-06,2.47E-06,1.02563,0.0131152,0.00637704,2.17E-06,2.14E-06,0.0246741,0.0168135,0.0235533,0.0130626,2.16E-06,0.0545531,0.00797507} 
    };

    //A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V
    private static final double[] MIXTURE_PROBS =
            new double[]{0.176513,0.0270202,0.0301127,0.0669246,0.207622,0.034662,0.0868259,0.0358616,0.03427,0.0428319,0.0466614,0.047875,0.0283695,0.0593123,0.0233828,0.0226822,0.00898452,0.00710292,0.00582299,0.00716226};
    
    private static final String IDENTIFIER = "R3.20C"; 
    
    public SAMRecode320CompPenalty(){
        super(MIXTURE_PROBS, ALPHAS, IDENTIFIER);
    } 

    public SAMRecode320CompPenalty(String justForTesting){// juit trying it out with Asi'fs numbers
        super(SAMMixures.RECODE3_20COMP_MIXTURE, SAMMixures.RECODE3_20COMP_FREQUENCIES, IDENTIFIER);
    } 



    public static void main(String[] args){
        System.out.println("testing sam mixtures");

        Fitness f = new Fitness(new double[]{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0});
        DirichletMixPenalty penalty = new SAMRecode320CompPenalty("just tryying this out");

        double p = penalty.calculate(f.get());
        System.out.println("penalty value "+p);
        
        System.out.println("End");
    }
}
