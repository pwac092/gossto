/*This file is part of GOssTo.
 GOssTo is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 GOssTo is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with GOssTo.  If not, see <http://www.gnu.org/licenses/>.
 */
package HSM;

import GOtree.Assignment;
import GOtree.GOTerm;
import Jama.Matrix;
import java.io.IOException;
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
    public Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.geneWiseSimilarityByMaximum(ontology);
    }

    @Override
    public Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        assert (ontology >= 0 && ontology < 3);

        float M = Float.NEGATIVE_INFINITY; //getting M

        final int N = numGOtermsPerOntology[ontology];
        Matrix result = new Matrix(N, N);

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                float dJiang = (float) -Math.log(lowestCommonAncestor(matrixAxis[ontology][i].getAncestors(), matrixAxis[ontology][j].getAncestors(), ontology));
                result.set(i, j, dJiang);
                result.set(j, i, dJiang);
                M = Math.max(M, dJiang);
            }
        }
        //Jiang normalisation:
        float[] normalizedDiagonal = new float[N];
        final float invM = 1.0f / M;
        for (int i = 0; i < N; i++) {
            normalizedDiagonal[i] = result.get(i, i) * invM;
        }

        float maxJiang = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                //set and calculate the Jiang value
                float value = (-2.0f * (result.get(i, j) * invM) + normalizedDiagonal[i] + normalizedDiagonal[j]);
                result.set(i, j, value);
                result.set(j, i, value);
                maxJiang = Math.max(maxJiang, value);
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                result.set(i, j, 1.0f - (result.get(i, j) / maxJiang)); //Jiang finishing normalisation
                result.set(j, i, result.get(i, j));
            }
        }

        logwriter.log("Completed HSM for " + shortOntologyName[ontology]);
        System.out.println("Completed Jiang for Ontology : " + longOntologyName[ontology]);
        return result;
    }
}
