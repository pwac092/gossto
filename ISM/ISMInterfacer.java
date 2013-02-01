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
package ISM;

import GOtree.Assignment;
import GOtree.GOTerm;
import ISM_ImplementationStrategies.ISM_validImplementation;
import java.io.IOException;
import java.util.ArrayList;
import org.apache.commons.math3.linear.RealMatrix;
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
    public RealMatrix getISMs(RealMatrix[] adjMatrices, GOTerm[][] matrixAxis, RealMatrix HSM, Assignment annotations,
            ArrayList<GOTerm> userProvidedTerms, String[] GO_relations, String dagChoice, int matrix, TinyLogger logger)
            throws IOException {

        ISM_validImplementation ism = new ISM_validImplementation(matrixAxis[matrix], HSM, GO_relations, annotations, true, false);
        RealMatrix result = ism.computeISM();

        return result;
    } 

    //get genewise ism data, the parameters for this method match those of the ISM Implementation constructor. For more detail look in the relevant class.
    public RealMatrix getGeneISMs(RealMatrix[] adjMatrices, GOTerm[][] matrixAxis, RealMatrix HSM, Assignment annotations,
            ArrayList<GOTerm> userProvidedTerms, String[] GO_relations, String dagChoice, int matrix, TinyLogger logger, boolean weightedJaccard)
            throws IOException {

        ISM_validImplementation ism = new ISM_validImplementation(matrixAxis[matrix], HSM, GO_relations, annotations, false, weightedJaccard);
        RealMatrix result = ism.computeISM();

        return result;
    }
}
