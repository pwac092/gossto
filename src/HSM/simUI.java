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
import HSM.GraphSimilarities.GraphSimilarity;
import HSM.GraphSimilarities.SimUISimilarity;
import Jama.Matrix;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron, Alfonso E. Romero
 *
 */
//Implements the simUI semantic similarity measure
public class simUI extends HSM {
    
    public simUI(GOTerm[] allTerms, String[] genes, String[][] goIds, GOTerm[][] axis, String[] targets, Assignment annotations, String[] relations, TinyLogger logw) {
        super(allTerms, genes, goIds, axis, targets, annotations, relations, logw);
        isAGraphBasedMeasure = true;    
    }

    @Override
    public Matrix calculateGeneWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        return super.calculateGraphGeneWiseSemanticSimilarity(ontology, new SimUISimilarity(annotations));
    }

    @Override
    public Matrix calculateTermWiseSemanticSimilarity(int ontology) throws IOException, OutOfMemoryError {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
