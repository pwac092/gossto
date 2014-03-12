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
//Implements Lin's semantic similarity measure
public class Lin extends HSM {

    public Lin(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
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
                float linTop = 0.0f - (float) Math.log(lowestCommonAncestor(matrixAxis[ontology][i].getAncestors(), matrixAxis[ontology][j].getAncestors(), ontology));
                result.set(i, j, linTop);
                result.set(j, i, linTop);

                if (linTop > M && linTop != 0) {
                    M = linTop; //get largest value for normalisation
                }
            }
        }

        //Lin normalisation:
        float [] normalizedDiagonal = new float [N];
        final float invM = 1.0f / M;
        for (int i = 0; i < N; i++) {
            normalizedDiagonal[i] = result.get(i, i) * invM;
        }

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                //create the bottom of the lin equation using the extracted values and a 0.001 to prevent a division by 0
                float linBottom = normalizedDiagonal[i] + normalizedDiagonal[j] + 0.001f;
                //create the top of the lin equation, divide it by the bottom and save the resultant lin value
                float val = 2.0f * (result.get(i, j) * invM) / linBottom;
                result.set(i, j, val);
                result.set(j, i, val);
            }
        }

        logwriter.log("Completed HSM for " + shortOntologyName[ontology]);
        System.out.println("Completed Lin for Ontology : Biological Process" + longOntologyName[ontology]);

        return result;
    }
}