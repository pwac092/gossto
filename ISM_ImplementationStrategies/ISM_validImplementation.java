/*This file is part of GOSSTO.
 GOSSTO is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GOSSTO is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GOSSTO.  If not, see <http://www.gnu.org/licenses/>.
 */
package ISM_ImplementationStrategies;

import GOtree.Assignment;
import GOtree.GOTerm;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

/**
 *
 * @author pwac092 This class computes the ISM for any termwise HSM provided. If
 * more than one tree is to be calculated, they have to be provided one at a
 * time.
 */
public class ISM_validImplementation {

    /*
     * P = 0; //transition probabilities
     * P(c,c) = 1 if node c is a leaf 
     * P(v,c) = (1 - N_v* / N_v) N_c/(Sum{u: v->u} N_u_)
     * W = I; //random walkers
     * epsilon = 0.001
     */
    /*utils*/
    //number of GOterms in the adjacency matrix
    private GOTerm[] subGoTerms;  //this holds all the goterms with annotations
    private String[] relations; //this holds the relations
    private HashMap<Integer, Integer> goTermIndex; //correspondence between goterms and their indices in the matrix.
    private Set<Integer> leafs; //just to cache the leafs.
    private int[] allIndices; //just to store them. We will definitely need them
    private int[] leafIndices; //again, we will definitely use them
    private HashMap<String, Integer> proteinIndices; //to store the protein indices for genewise calculations
    private double maxNumberOfAnnotations; //an integer to store the maximnun number of annotations of any node in the tree
    /*ISM elements*/
    private OpenMapRealMatrix P; //size: V
    private OpenMapRealMatrix W;
    private RealMatrix RWC;
    private RealMatrix ISM;
    private double epsilon;
    /*HSM*/
    private RealMatrix HSM;
    /*Annotations*/
    private Assignment annotations; //the actual annotations. 
    private HashMap<Integer, Integer> numAnnotations; //just a cache for regular annotations
    private HashMap<Integer, Integer> numAnnotationsStar; //just a cache for Star (i.e not in children) values
    /*We will produce both genewise and termwise from this class. 
     * this will reduce the already massive code footprint.
     */
    private boolean termwise;
    private boolean weightedJaccard;

    public ISM_validImplementation(GOTerm[] ISM_currentGoTerms, RealMatrix HSM, String[] ISM_currentRelations, Assignment ISM_Annotations, boolean termwise, boolean wJaccard) {


        //0. Utils

        //0.1 various caches for speedup
        this.annotations = ISM_Annotations;//copy the annotations
        this.numAnnotations = new HashMap<Integer, Integer>(); //fire up the cache
        this.numAnnotationsStar = new HashMap<Integer, Integer>(); //fire up the cache      

        //0.2 load the affected goterms
        this.subGoTerms = ISM_currentGoTerms;
        
        //0.3 just some variables to keep useful data.
        this.relations = ISM_currentRelations;
        this.maxNumberOfAnnotations = this.getMaxOntology();

        //0.4 HSM
        this.HSM = HSM.copy();

 
        //0.5.0 fire up the cache for the indices and load it up
        this.goTermIndex = new HashMap<Integer, Integer>();
        //now we load all the indices for the goterms.
        for (int i = 0; i < this.getNumGoTerms(); i++) {
            this.goTermIndex.put(this.subGoTerms[i].getNumericId(), i);
        }
        //0.5.1 now that the goterms are set, we can load the indices.
        this.allIndices = this.getAllIndices();


        //0.6 finding the proper leafs.
        this.leafs = new HashSet<Integer>(); //all leafs are here
        this.setAllLeafs();
        this.leafIndices = this.getLeafIndices();

        //0.7 ISM
        //to differentiate between genewise and termwise.
        this.termwise = termwise;
        //to differentiate between weighted and unweighted jaccard
        this.weightedJaccard = wJaccard;


        //we load the indices for the annotations only in the case
        //were we compute the genewise similarity
        if (!this.termwise) {
            //this is just to index all the annotations.
            this.proteinIndices = new HashMap<String, Integer>();
            for (String uniqueAnnotation : this.annotations.getRowIdentifiers()) {
                int proteinIndex;
                //the protein already exists in the list.
                if (!proteinIndices.containsKey(uniqueAnnotation)) {
                    proteinIndex = this.proteinIndices.size();
                    this.proteinIndices.put(uniqueAnnotation, proteinIndex);
                }
            }
        }


        //0.8 RWC
        //Different sets have to be traversed
        if (this.termwise) {
            this.RWC = new OpenMapRealMatrix(this.getNumGoTerms(), this.getNumGoTerms());
        } else {
            this.RWC = new OpenMapRealMatrix(this.proteinIndices.size(), this.proteinIndices.size());
        }
        this.ISM = null; //new OpenMapRealMatrix(this.getNumGoTerms(), this.getNumGoTerms());


        //1. initialise transitionprobabilities
        //we use a sparse matrix, so we don't need to put zeroes anywhere.
        this.P = new OpenMapRealMatrix(this.getNumGoTerms(), this.getNumGoTerms());
        //2. initialise random walkers.
        //keep in mind that this makes sense because of the way the indexes were
        //loaded into this.goTermIndex. Otherwise, we would have to retrieve 
        //from this.goTermIndex
        this.W = new OpenMapRealMatrix(this.getNumGoTerms(), this.getNumGoTerms());
        for (int i = 0; i < this.getNumGoTerms(); i++) {
            this.W.setEntry(i, i, 1.0);
        }
        //3.
        //just set the convergence limit.
        this.epsilon = 0.001;

    }

    public RealMatrix computeISM() {
        //Step 0. Initialise transition probabilities
        this.initialiseTransitionProbabilities();
        //Step 1. Walk!
        this.walk();
        //Step 2. Compute the random wal contribution
        if (true == this.termwise) {
            this.setRandomWalkContributionTermwise();
        } else {
            this.setRandomWalkContributionGenewise();
        }
        //Step 3. get the ISM which we return.
        this.setISM();
        return this.ISM;
    }

    private void initialiseTransitionProbabilities() {

        //0.1 check if the node is a leaf, if it is, put a 1 into it.
        for (GOTerm currentGoTerm : this.subGoTerms) {
            if (this.leafs.contains(currentGoTerm.getNumericId())) {
                int leafIndex = this.goTermIndex.get(currentGoTerm.getNumericId());
                this.P.setEntry(leafIndex, leafIndex, 1.0);
                continue;
            }
        }
        for (GOTerm currentGoTerm : this.subGoTerms) {
            int N_v = this.getNumberOfAnnotations(currentGoTerm);
            int N_vStar = this.getNumberOfAnnotationsStar(currentGoTerm);

            this.setTransitionProbabilitiesNonLeaf(currentGoTerm, N_v, N_vStar);
        }
    }

    private void setTransitionProbabilitiesNonLeaf(GOTerm currentGoTerm, int N_v, int N_vStar) {
        //P(v,c) = (1 - N_v* / N_v) N_c/(Sum{u: v->u} N_u_)
        //P(v,c) = A * N_C/B

        if (N_v == N_vStar) {
            // The fact that N_v is equal to N_vStar 
            // gives a 0.0 in the matrix entries, because of the leading factor
            // 1-N_v/N_vStar, which does not to be specified due to the sparse
            // matrix used...
            return;
        }

        double A = (1.0 - ((double) N_vStar / (double) N_v));

        //we get all the children
        //now, this is a tricky one.
        //the way I understand it so far is that all relations are the same, that is, 
        //the random walker will not differentiate between a part_of and a is_a relation
        //so, we get ALL the children, adding them all up in a list, disregarding the relations
        //completely. This is viable, since the random walker would be able to reach that 
        //node following either one of the relations.
        //children will the contain ALL the children, following every possible relation

        //the number of children per node is not very large, ArrayList is O(1) insertion, and O(n)
        //so it's ok.
        Set<GOTerm> children = new HashSet<GOTerm>();
        for (String currentRelation : this.relations) {
            children.addAll(currentGoTerm.getChildrenForRelation(currentRelation));
        }

        //we count all the annotations in the children
        //B keesp the sum of N_u for every child.
        //this number is the same given eery child of v.
        int N_u = 0;
        for (GOTerm currentChild : children) {
            N_u += this.getNumberOfAnnotations(currentChild);
        }
        //we need to use two loops. First we compute the total sum of 
        //annotations in the children, and the we modify the matrix P.
        //keep in mind that these nodes are allways non leaf nodes, therefore, B 
        //should not remain zero unless the children are actually not annotating
        //any genes.
        //this is why we specify an initial value for B. This will put 
        //a very low transition probability to that node.

        //P(v,c) = A * N_C/B

        final double inv_N_u = 1.0 / N_u;

        Map<Integer, GOTerm> sortedChildren = new TreeMap<Integer, GOTerm>();

        for (GOTerm currentChild : children) {
            sortedChildren.put(currentChild.getNumericId(), currentChild);
        }


        for (int currentChildId : sortedChildren.keySet()) {
            GOTerm currentChild = sortedChildren.get(currentChildId);

            int N_c = this.getNumberOfAnnotations(currentChild);

            //we check whether the annotations are in excess of zero. 
            //if the node does not have any annotations, then there is no need 
            //in putting anything in it, since the structure is a sparse matrix.


            if (N_c > 0) {
                double newEntry = A * (N_c * inv_N_u);
                int v = this.goTermIndex.get(currentGoTerm.getNumericId());
                int c = this.goTermIndex.get(currentChild.getNumericId());
                this.P.setEntry(c, v, newEntry);
            }
        }
    }

    private void walk() {

        OpenMapRealMatrix W_star = new OpenMapRealMatrix(this.W);
        double convergence;
        do {
            this.W = W_star.copy();
            W_star = this.P.multiply(this.W).copy();
            convergence = W_star.subtract(this.W).getFrobeniusNorm();
        } while (convergence > this.epsilon);
        this.W = W_star.copy();
    }

    private void setRandomWalkContributionTermwise() {

        System.out.println("Submatrix (RWC)");
        this.RWC = this.W.getSubMatrix(this.leafIndices, this.allIndices);
        System.out.println("Transpose (RWC)");

        this.RWC = this.RWC.transpose();
        System.out.println("Submatrix (HSM)");

        RealMatrix subHSM = this.HSM.getSubMatrix(this.leafIndices, this.leafIndices);
        System.out.println("RWC * HSM");

        this.RWC = this.RWC.multiply(subHSM);
        System.out.println("Submatrix (W)");

        RealMatrix subW = this.W.getSubMatrix(this.leafIndices, this.allIndices);
        System.out.println("RWC * SubMatrixW");

        this.RWC = this.RWC.multiply(subW);
    }

    private void setRandomWalkContributionGenewise() {
        //0. get matrix A
        RealMatrix A = this.getMatrixA();
        //1. multiply both matrices.
        RealMatrix B = this.W.getSubMatrix(this.leafIndices, this.allIndices).multiply(A);
        //2. calculate the RWC
        //2.0 traverse all the products.
        //set the value for all eht rows for this column and this row
        //for RWC column_index  == row_index
        for (int i = 0; i < this.RWC.getRowDimension(); i++) {
            for (int j = i; j < this.RWC.getColumnDimension(); j++) {
                double jaccardIndex = this.getJaccardIndex(B.getColumn(i), B.getColumn(j));
                this.RWC.setEntry(i, j, jaccardIndex);
                this.RWC.setEntry(j, i, jaccardIndex);
            }
        }
    }

    private RealMatrix getMatrixA() {

        RealMatrix A = new OpenMapRealMatrix(this.getNumGoTerms(), this.annotations.sizeGenes());
        for (GOTerm currentGoTerm : this.subGoTerms) {
            //0. check for NStar value > 0, since this indicates there
            //is an annotation 
            if (this.getNumberOfAnnotationsStar(currentGoTerm) > 0) {
                //0. Get all the genes annotating the current node
                //this will hodl the difference of the parent set with the child set
                Set<String> uniqueAnnotations = new HashSet<String>(this.annotations.getProteinsForGOTerm(currentGoTerm.getGOid()));
                //1. Get all the genes annotating the children of the current node
                //1.0 get all the children
                Set<GOTerm> children = new HashSet<GOTerm>();
                for (String currentRelation : this.relations) {
                    children.addAll(currentGoTerm.getChildrenForRelation(currentRelation));
                }
                //1.1 get all the genes for the children.
                Set<String> childrenAnnotations = new HashSet<String>();
                for (GOTerm currentChild : children) {
                    childrenAnnotations.addAll(this.annotations.getProteinsForGOTerm(currentChild.getGOid()));
                }
                //2. Traverse the difference and count the number of terms annotating this et of genes.
                //2.0 obtain the difference between the two nodes. that is, the annotations
                //that are unique to the current node..
                uniqueAnnotations.removeAll(childrenAnnotations);
                for (String uniqueAnnotation : uniqueAnnotations) {
                    //this is a tricky one. We are not sure what "directly" means in the paper.
                    //but it should not be a very complicated problem to solve.
                    int count = this.annotations.getGOTermScoresForProteinId(uniqueAnnotation).keySet().size();
                    //0. get the protein id.
                    A.setEntry(this.goTermIndex.get(currentGoTerm.getNumericId()), this.proteinIndices.get(uniqueAnnotation), 1.0 / count);
                }
            }
        }
        return A;
    }

    private double getJaccardIndex(double[] distributionProductOne, double[] distributionProductTwo) {
        //0. gather data
        double combinedSum = 0;
        double sumProductOne = 0, sumProductTwo = 0;
        double IC = 1.0;
        for (int i = 0; i < this.leafIndices.length; i++) {

            //we need to fetch  the information content of the nodes if we use weighted jaccard.
            if (this.weightedJaccard) {
                IC = -Math.log(this.numAnnotations.get(this.leafIndices[i]) / this.maxNumberOfAnnotations);
            }
            //compute data independently
            combinedSum += distributionProductOne[i] * distributionProductTwo[i] * IC;
            sumProductOne += distributionProductOne[i] * IC;
            sumProductTwo += distributionProductTwo[i] * IC;
        }
        //1. compute the jaccard index
        return (combinedSum / (sumProductOne + sumProductTwo - combinedSum));
    }

    private void setISM() {
        this.ISM = (this.RWC.add(HSM)).scalarMultiply(0.5);
    }

    //**************************************************************************
    //**************************************************************************
    //UTILITY FUNCTIONS
    //**************************************************************************
    //**************************************************************************
    private int getNumberOfAnnotations(GOTerm currentGoTerm) {
        //look for the annotations in the cache, if we find it, return it.
        int id = currentGoTerm.getNumericId();
        if (this.numAnnotations.containsKey(id)) {
            return this.numAnnotations.get(id);
        } else { //this number is not cached yet.
            int num = this.annotations.countNumberOfGenesForGOTerm(currentGoTerm.getGOid());
            this.numAnnotations.put(id, num);
            return num;
        }
    }

    private int getNumberOfAnnotationsStar(GOTerm currentGoTerm)
            throws IllegalArgumentException {

        //look for the annotations in the cache, if we find it, return it.
        int id = currentGoTerm.getNumericId();
        if (this.numAnnotationsStar.containsKey(id)) {
            return this.numAnnotationsStar.get(id);
        } else { //this number is not cached yet, so we need to compute it.

            //0. get number of annotations for this node.
            int currentGoTermAnnotationCount = this.getNumberOfAnnotations(currentGoTerm);
            //1. get all annotations for the children.

            Set<GOTerm> children = new HashSet<GOTerm>();
            for (String currentRelation : this.relations) {//again, we consider all relations at once.
                children.addAll(currentGoTerm.getChildrenForRelation(currentRelation));
            }


            Set<String> proteinsAnnotatedToChildren = new HashSet<String>();
            for (GOTerm child : children) {
                proteinsAnnotatedToChildren.addAll(this.annotations.getProteinsForGOTerm(child.getGOid()));
            }

            //there is no need to actually make the set difference. Once we 
            //get the number of annotations in the parent, we should just
            //substract to that number the number of unique annotations in the
            //children. This would be the number of annotations in the parent
            //that belong to none of the children.

            int childrenAnnotationCount = proteinsAnnotatedToChildren.size();

            int retVal = currentGoTermAnnotationCount - childrenAnnotationCount;

            this.numAnnotationsStar.put(id, retVal);

            return retVal;
        }
    }

    private void setAllLeafs() {

        for (GOTerm term : this.subGoTerms) {
            boolean isLeaf = true;
            for (String relation : this.relations) {
                //is it a leaf for ALL the ontologies? otherwise, 
                //it is not really a leaf..
                if (!this.isLeafISM(term, relation)) {
                    isLeaf = false;
                    break;
                }
            }
            if (isLeaf) {
                this.leafs.add(term.getNumericId());
            }
        }
    }

    /*we write this to conform to the defintion of leaf in ISM..
     * Basically, we need to find if the children of the current node belongs 
     * to the "ontology". By "ontology" we mean the subset of nodes with annotations
     * The easies way of doing this is: fetch the children, iterate through them
     * and find the number of annotations. If it all of the children do not have annotations
     * then it is a leaf, oftherwise, it is not.
     */
    private boolean isLeafISM(GOTerm currentTerm, String relation) {

        List<GOTerm> children = currentTerm.getChildrenForRelation(relation);

        //if it is emtpy, then the thing is a leaf, for sure
        //keep in mind that we only check children for nodes with annotations
        if (children.isEmpty()) {
            return true;
        }
        //check the annotation count for each of the children, in case the thing
        //has children

        for (GOTerm currentChild : children) {
            //the easy way: is the intersection between children and this.goTerms 
            //emtpy?
            if (this.getNumberOfAnnotations(currentChild) > 0) {
                return false;
            }
        }
        return true;
    }

    private int[] getLeafIndices() {
        int allIdx[] = new int[this.leafs.size()];
        int j = 0;

        for (int leafNumericId : this.leafs) {
            allIdx[j++] = this.goTermIndex.get(leafNumericId);
        }
        Arrays.sort(allIdx); // could this be a problem?
        return allIdx;
    }

    private int[] getAllIndices() {
        int[] allIdx = new int[this.goTermIndex.size()];
        for (int i = 0; i < this.goTermIndex.size(); i++) {
            allIdx[i] = i;
        }
        return allIdx;
    }

    //Getters yeah..this doesn't work here...
    //inline int getNumGoTerms(){return this.numGoTerms};
    private int getNumGoTerms() {
        return this.subGoTerms.length;
    }
    
    private double getMaxOntology()
   {
       double max = Double.NEGATIVE_INFINITY;
       for (GOTerm term : this.subGoTerms) {
           max = Math.max(max, this.annotations.countNumberOfGenesForGOTerm(term.getGOid()));
       }
       return max;
   }
 
    
}
