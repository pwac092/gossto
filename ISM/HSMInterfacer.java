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

/**
 *
 * @author Samuel Heron
 */
import GOtree.GOTerm;
import HSM.HSM;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Set;
import org.apache.commons.math3.linear.OpenMapRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import util.TinyLogger;

/**
 * This class retreives and interfaces with HSM instances, acting like a buffer
 * between their classes and the main ISM.java class
 */
public class HSMInterfacer {

    /**
     * Used for writing messages to the log file
     */
    private TinyLogger logwriter;
    /**
     * Stores retrieved HSM instance
     */
    private HSM chosenHSM;
    private final Set<GOTerm> targets;
    private final GOTerm[][] matrixAxis;

    /**
     * Constructor: Instantiates the log file variable if a log file is to be
     * written
     */
    public HSMInterfacer(TinyLogger logw, Set<GOTerm> targets, GOTerm[][] matrixAxis) {
        this.logwriter = logw;
        this.targets = targets;
        this.matrixAxis = matrixAxis;
    }

    public String[] getComputedGenes() {
        return this.chosenHSM.getSubSetGenes();
    }

    public boolean isAGraphBasedMeasure() {
        return chosenHSM.isAGraphBasedMeasure();        
    }

    //Retrieves the specific HSm instance, parameters detailed above the 'getHSMinstance()' method
    public void retrieveHSMinstance(String name, Object[] params) throws IOException {
        HSM hsmInstance = null;
        // we search for the class which name is the String 'name'
        Class<? extends HSM> MyClass;
        try {
            MyClass = Class.forName("HSM." + name).asSubclass(HSM.class);
            // fetches all relevant constructors for the specified class

            Constructor<?>[] constructors = MyClass.getConstructors();
            hsmInstance = (HSM) constructors[0].newInstance(params);

        } catch (InstantiationException ex) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " is not instantiable.");
            System.err.println("ERROR: class " + name + " is not instantiable.");
            System.exit(-1);
        } catch (IllegalAccessException ex) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " has an inaccessible constructor.");
            System.err.println("ERROR: class " + name + " has an inaccessible constructor.");
            System.exit(-1);
        } catch (IllegalArgumentException e) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " has an illegal argument.");
            System.err.println("ERROR: class " + name + " has an illegal argument.");
            System.exit(-1);
        } catch (InvocationTargetException e) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " could not be invoked.");
            e.getCause().printStackTrace(System.err);
            System.err.println("ERROR: class " + name + " could not be invoked.");
            System.exit(-1);
        } catch (ClassNotFoundException e) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " not found.");
            System.err.println("ERROR: class " + name + " not found.");
            System.exit(-1);
        } catch (SecurityException e) {
            this.logwriter.logAndCloseWriter("############ERROR: class " + name + " has triggered a security exception.");
            System.err.println("ERROR: class " + name + " has triggered a security exception.");
            System.exit(-1);
        }
        this.chosenHSM = hsmInstance;
    }

    public RealMatrix returnGeneWiseResults(int matrix) throws IOException {
        return this.chosenHSM.calculateGeneWiseSemanticSimilarity(matrix);
    }

    //Retrieves the HSM results, the parameter specifying whether we want to force it to return the gene simiarity results (only required fro printing)
    public RealMatrix returnTermWiseResults(int matrix) throws IOException {
        if (this.targets == null || this.targets.isEmpty()) {
            return chosenHSM.calculateTermWiseSemanticSimilarity(matrix);
        } else {
            return returnTrimmedMatrix(this.chosenHSM.calculateTermWiseSemanticSimilarity(matrix), matrix);
        }
    }

    private RealMatrix returnTrimmedMatrix(RealMatrix in, int matrix) {
        int size = 0, rowInd = 0, colInd = 0;
        RealMatrix trimmedMatrix;
        trimmedMatrix = null;

        for (int i = 0; i < in.getRowDimension(); i++) {
            if (targets.contains(matrixAxis[matrix][i])) {
                size++;
            }
        }
        if (size > 0) {
            trimmedMatrix = new OpenMapRealMatrix(size, size);
            for (int i = 0; i < in.getRowDimension(); i++) {
                if (targets.contains(matrixAxis[matrix][i])) {
                    for (int j = 0; j < in.getRowDimension(); j++) {
                        if (targets.contains(matrixAxis[matrix][j])) {
                            trimmedMatrix.setEntry(rowInd, colInd, in.getEntry(i, j));
                            colInd++;
                        }
                    }
                    colInd = 0;
                    rowInd++;
                }
            }

        }
        return trimmedMatrix;
    }
}
