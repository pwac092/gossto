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
package HSM;

import GOtree.Assignment;
import GOtree.GOTerm;
import java.io.IOException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//Implements Jiang & Conrath's semantic similarity measure
public class Jiang extends HSM {

    public Jiang(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
    }

    @Override
    public RealMatrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.geneWiseSimilarityByMaximum(ontology);
    }

    @Override
    public RealMatrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        assert (ontology >= 0 && ontology < 3);

        double M = Double.NEGATIVE_INFINITY; //getting M

        final int N = numGOtermsPerOntology[ontology];
        RealMatrix result = new Array2DRowRealMatrix(N, N);

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                double dJiang = -Math.log(lowestCommonAncestor(matrixAxis[ontology][i].getAncestors(), matrixAxis[ontology][j].getAncestors(), ontology));
                result.setEntry(i, j, dJiang);
                result.setEntry(j, i, dJiang);
                M = Math.max(M, dJiang);
            }
        }
        //Jiang normalisation:
        double[] normalizedDiagonal = new double[N];
        final double invM = 1.0 / M;
        for (int i = 0; i < N; i++) {
            normalizedDiagonal[i] = result.getEntry(i, i) * invM;
        }

        double maxJiang = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                //set and calculate the Jiang value
                double value = (-2.0 * (result.getEntry(i, j) * invM) + normalizedDiagonal[i] + normalizedDiagonal[j]);
                result.setEntry(i, j, value);
                result.setEntry(j, i, value);
                maxJiang = Math.max(maxJiang, value);
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                result.setEntry(i, j, 1.0 - (result.getEntry(i, j) / maxJiang)); //Jiang finishing normalisation
                result.setEntry(j, i, result.getEntry(i, j));
            }
        }

        logwriter.log("Completed HSM for " + shortOntologyName[ontology]);
        System.out.println("Completed Jiang for Ontology : " + longOntologyName[ontology]);
        return result;
    }
}
