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
//Implements Resnik's semantic similarity measure
public class Resnik extends HSM {
    
    public Resnik(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
    }

    @Override
    public Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        assert (ontology >= 0 && ontology < 3);
        final int N = numGOtermsPerOntology[ontology];
        
        super.logwriter.showMessage("Size of the semsim matrix: " + N);
        super.logwriter.showMessage("Size of the semsim matrix (in MBs): " + N*N*4.0/(1024*1024));
        
        Matrix result = new Matrix(N, N);

        for (int i = 0; i < N; i++) {
            for (int j = i; j < N; j++) {
                float res = (float) -Math.log(lowestCommonAncestor(matrixAxis[ontology][i].getAncestors(), matrixAxis[ontology][j].getAncestors(), ontology));
                result.set(i, j, res);
                result.set(j, i, res);
            }
        }

        logwriter.showTimedMessage("Completed Resnik for Ontology " + longOntologyName[ontology]);

        return result;
    }

    @Override
    public Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.geneWiseSimilarityByMaximum(ontology);
    }
}
