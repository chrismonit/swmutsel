package swmutsel.model;

import cern.colt.Arrays;
import org.apache.commons.math3.special.Gamma;
import swmutsel.utils.CoreUtils;

/**
 *
 * @author Christopher Monit c.monit.12@ucl.ac.uk
 */
public abstract class DirichletMixPenalty implements Penalty {
    protected double[] MIXTURE_PROBS;
    protected double[][] ALPHAS;
    protected String IDENTIFIER;
    
    protected double[] invBetaTerms;
    
    public DirichletMixPenalty(double[] mix_probs, double[][] alphas, String identifier){
        this.MIXTURE_PROBS = mix_probs;
        this.ALPHAS = alphas;
        this.IDENTIFIER = identifier;
        
        this.invBetaTerms = calculateInvBetaTerms();
    }

    
    private final double[] calculateInvBetaTerms(){ // return { 1/B({alphas}) }
        double[] terms = new double[MIXTURE_PROBS.length];
        for (int iMix = 0; iMix < MIXTURE_PROBS.length; iMix++) {
            double sumAlpha = CoreUtils.sum(this.ALPHAS[iMix]);

            double productGammaAlpha = 1.0;
            for (int iAlpha = 0; iAlpha < ALPHAS[0].length; iAlpha++)
                productGammaAlpha *= Gamma.gamma(this.ALPHAS[iMix][iAlpha]);
            
            terms[iMix] = Gamma.gamma(sumAlpha)/productGammaAlpha;
        }// iMix

        return terms;
    }
    

    @Override
    public final double calculate(final double[] parameters){
        double[] frequencies = calculateFrequencies(parameters); 
        
        double sum = 0.0; // sum over mixture components
        System.out.println(Arrays.toString(invBetaTerms));
        for (int iMix = 0; iMix < this.MIXTURE_PROBS.length; iMix++) {
            double prod = 1.0;
            for (int iFreq = 0; iFreq < frequencies.length; iFreq++)
                prod *= Math.pow(frequencies[iFreq], this.ALPHAS[iMix][iFreq]-1. );
            System.out.printf("mixprob %f, invBeta %f, prod %f %n", this.MIXTURE_PROBS[iMix], this.invBetaTerms[iMix], prod);
            sum += this.MIXTURE_PROBS[iMix] * this.invBetaTerms[iMix] * prod;
        }// iMix
        System.out.println("caluclate, sum"+sum);
        return Math.log(sum); 
    }
    

    
    // TODO this is the same as SwMutSel set/get getCodonFrequencies - could reuse that?
    public static final double[] calculateFrequencies(double[] fitness){
        double[] frequencies = new double[fitness.length];
        double z = 0.0;
        for (int i = 0; i < fitness.length; i++) {
            frequencies[i] = Math.exp(fitness[i]);
            z += frequencies[i];
        }

        for (int i = 0; i < fitness.length; i++) frequencies[i] /= z;

        return frequencies;
    }

}
