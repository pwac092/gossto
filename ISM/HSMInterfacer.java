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

/**
 *
 * @author Samuel Heron
 */
import GOtree.GOTerm;
import HSM.HSM;
import Jama.Matrix;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import util.TinyLogger;

/**
 * This class retreives and interfaces with HSM instances, acting like a buffer
 * between their classes and the main ISM.java class
 */
public class HSMInterfacer {

    /**
     * Used for writing messages to the log file
     */
    private final TinyLogger logwriter;
    /**
     * Stores retrieved HSM instance
     */
    private HSM chosenHSM;
    private final Set<GOTerm> targets;
    private final GOTerm[][] matrixAxis;
    private final String[] targetGenes;
    private Matrix originalMatrix;

    /**
     * Constructor: Instantiates the log file variable if a log file is to be
     * written
     */
    public HSMInterfacer(TinyLogger logw, Set<GOTerm> targets, GOTerm[][] matrixAxis, String[] targetGenes) {
        this.logwriter = logw;
        this.targets = targets;
        this.matrixAxis = matrixAxis;
        this.targetGenes = targetGenes;
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

    public Matrix getOriginalCachedMatrix() {
        return this.originalMatrix;
    }

    public Matrix returnGeneWiseResults(int matrix) throws IOException {
        if (chosenHSM.getNumGOTermsPerOntology(matrix) == 0) {
            // this case might happen if the organism has no annotation in that ontology
            return null;
        }

        this.originalMatrix = this.chosenHSM.calculateGeneWiseSemanticSimilarity(matrix);
        if (this.targetGenes == null || this.targetGenes.length == 0) {

            return this.originalMatrix;
        } else {
            return returnTrimmedMatrixForGenes(this.originalMatrix);
        }
    }

    //Retrieves the HSM results, the parameter specifying whether we want to force it to return the gene simiarity results (only required fro printing)
    public Matrix returnTermWiseResults(int matrix) throws IOException {
        if (chosenHSM.getNumGOTermsPerOntology(matrix) == 0) {
            // this case might happen if the organism has no annotation in that ontology
            return null;
        }
        
        this.originalMatrix = chosenHSM.calculateTermWiseSemanticSimilarity(matrix);

        if (this.targets == null || this.targets.isEmpty()) {
            return this.originalMatrix;
        } else {
            return returnTrimmedMatrix(this.originalMatrix, matrix);
        }
    }

    private Matrix returnTrimmedMatrixForGenes(Matrix in) {
        Matrix trimmedMatrix;
        trimmedMatrix = null;

        int size = 0, rowInd = 0, colInd = 0;

        Set<String> selGenez = new HashSet<String>();
        selGenez.addAll(Arrays.asList(this.targetGenes));
        String[] allGenes = this.chosenHSM.getSubSetGenes();

        for (int i = 0; i < in.getRowDimension(); i++) {
            if (selGenez.contains(allGenes[i])) {
                ++size;
            }
        }

        if (size > 0) {
            trimmedMatrix = new Matrix(size, size);
            for (int i = 0; i < in.getRowDimension(); i++) {

                if (selGenez.contains(allGenes[i])) {
                    for (int j = 0; j < in.getRowDimension(); j++) {
                        if (selGenez.contains(allGenes[j])) {
                            trimmedMatrix.set(rowInd, colInd, in.get(i, j));
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

    private Matrix returnTrimmedMatrix(Matrix in, int matrix) {
        int size = 0, rowInd = 0, colInd = 0;
        Matrix trimmedMatrix;
        trimmedMatrix = null;

        for (int i = 0; i < in.getRowDimension(); i++) {
            if (targets.contains(matrixAxis[matrix][i])) {
                size++;
            }
        }
        if (size > 0) {
            trimmedMatrix = new Matrix(size, size);
            for (int i = 0; i < in.getRowDimension(); i++) {
                if (targets.contains(matrixAxis[matrix][i])) {
                    for (int j = 0; j < in.getRowDimension(); j++) {
                        if (targets.contains(matrixAxis[matrix][j])) {
                            trimmedMatrix.set(rowInd, colInd, in.get(i, j));
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
