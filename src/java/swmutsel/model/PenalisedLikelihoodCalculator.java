package swmutsel.model;

import java.util.Arrays;
import pal.tree.Tree;
import swmutsel.model.parameters.Fitness;
import swmutsel.model.parameters.Parameter;
import swmutsel.utils.CoreUtils;
import swmutsel.utils.Pair;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Wraps LikelihoodCalculator and adds penalty on Fitness parameters to log-likelihood.
 *
 *
 * Author: Asif Tamuri (tamuri@ebi.ac.uk)
 * Date: 18/10/2013 21:22
 * 
 * Modified: Christopher Monit (c.monit.12@ucl.ac.uk)
 * Date: 7 Sep 2017 17:18:42 BST
 */
public class PenalisedLikelihoodCalculator extends LikelihoodCalculator {
    private final List<? extends Parameter> penalisedParameters;
    private final Penalty penalty;

    // Applying same penalty to all sets of fitnesses
    public PenalisedLikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, LinkedHashMap<String, SubstitutionModel> models, Penalty penalty, List<? extends Parameter> penalisedParameters) {
        super(tree, siteStates, models);
        this.penalisedParameters = penalisedParameters;
        this.penalty = penalty;
    }


    public PenalisedLikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, LinkedHashMap<String, SubstitutionModel> models, Penalty penalty, Parameter penalisedParameter) {
        this(tree, siteStates, models, penalty, Arrays.asList(penalisedParameter));
    }
    
    @SuppressWarnings("unchecked")
    public PenalisedLikelihoodCalculator(Tree tree, Map<String, Byte> siteStates, SubstitutionModel model, Penalty penalty, Parameter penalisedParameter) {
        this(tree, siteStates, CoreUtils.getLinkedHashMap(Pair.of("ALL", model)), penalty, penalisedParameter);
    }

    @Override
    public double getLogLikelihood() {
        return super.getLogLikelihood() + calculatePenalty();
    }

    @Override
    public double updateTree(Tree tree, Map<Integer, Integer> newToOldNumbering, Set<Integer> updateNodes) {
        return super.updateTree(tree, newToOldNumbering, updateNodes) + calculatePenalty();
    }

    @Override
    public double getLogLikelihoodForBranch(int node, double branchLength) {
        double out = super.getLogLikelihoodForBranch(node, branchLength);
        out += calculatePenalty();
        return out;
    }
    
    // NB the penalty.calculate returns logged value
    private double calculatePenalty() {
        double penalty = 0.0;
        for (Parameter param : this.penalisedParameters){
            if (param instanceof Fitness) {
                penalty += this.penalty.calculate(param.getOptimisable());
            }
        }
        return penalty;
    }
}
