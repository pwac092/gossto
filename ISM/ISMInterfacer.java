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
package ISM;

import GOtree.Assignment;
import GOtree.GOTerm;
import ISM_ImplementationStrategies.ISM_validImplementation;
import Jama.Matrix;
import java.io.IOException;
import java.util.ArrayList;
import util.TinyLogger;

/**
 *
 * @author Samuel Heron
 */
//This class interfaces with the ISM implementation classes & thus acts like a buffer
public class ISMInterfacer {

    //Nothing requires instantiation
    ISMInterfacer() {
    }

    //get termwise ism data, the parameters for this method match those of the ISM Implementation constructor. For more detail look in the relevant class.
    public Matrix getISMs(GOTerm[][] matrixAxis, Matrix HSM, Assignment annotations,
            ArrayList<GOTerm> userProvidedTerms, String[] GO_relations, String dagChoice, int matrix, TinyLogger logger)
            throws IOException {

        if (HSM == null) {
            return null;
        } else {

            ISM_validImplementation ism = new ISM_validImplementation(matrixAxis[matrix], HSM, GO_relations, annotations, true, false, logger);
            Matrix result = ism.computeISM();

            return result;
        }
    }

    //get genewise ism data, the parameters for this method match those of the ISM Implementation constructor. For more detail look in the relevant class.
    public Matrix getGeneISMs(GOTerm[][] matrixAxis, Matrix HSM, Assignment annotations,
            ArrayList<GOTerm> userProvidedTerms, String[] GO_relations, String dagChoice, int matrix, TinyLogger logger, boolean weightedJaccard)
            throws IOException {
        if (HSM == null) {
            return null;
        } else {
            ISM_validImplementation ism = new ISM_validImplementation(matrixAxis[matrix], HSM, GO_relations, annotations, false, weightedJaccard, logger);
            Matrix result = ism.computeISM();

            return result;
        }
    }
}
