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
package GOtree;
 
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
//import Jama.Matrix;
import java.io.*;
import java.util.*;

public class Assignment {

    /**
     * Integer identifier corresponding to each protein/gene ID
     */
    Map<String, Integer> indexByGene;
    /**
     * Integer identifier corresponding to each GO Term
     */
    Map<String, Integer> indexByGoTerm;
    /**
     *
     */
    Map<Integer, Integer> countByGoTerm;
    /**
     * String ID of a protein/gene corresponding to a integer identifier
     */
    List<String> geneById;
    /**
     * String ID of a GO term corresponding to a integer identifier
     */
    List<String> goTermById;
    /**
     * For each protein (identifier by its integer ID) and GO term, represents,
     * in a form of sparse matrix, the strength of the score of the association
     * using a real number
     */
    Map<Integer, Map<Integer, Double>> values;

    public Assignment() {
        this.indexByGene = new HashMap<String, Integer>();
        this.indexByGoTerm = new HashMap<String, Integer>();
        this.countByGoTerm = new HashMap<Integer, Integer>();
        this.geneById = new ArrayList<String>();
        this.goTermById = new ArrayList<String>();
        this.values = new HashMap<Integer, Map<Integer, Double>>();
    }

    public double getMax() {
        double maximum = Double.NEGATIVE_INFINITY;

        for (int key : this.values.keySet()) {
            for (int key2 : this.values.get(key).keySet()) {
                maximum = Math.max(maximum, this.values.get(key).get(key2));
            }
        }
        return maximum;
    }

    public static Assignment linearCombine(final Assignment a1, double alpha1, final Assignment a2, double alpha2) {
        Assignment result = a1.clone().multiplyValuesByScalar(alpha1);
        result.add(a2.multiplyValuesByScalar(alpha2));
        return result;
    }

    public boolean containsProtein(String proteinId) {
        return indexByGene.containsKey(proteinId);
    }

    public Assignment(final Assignment other) {
        this.indexByGene = new HashMap<String, Integer>(other.indexByGene);
        this.indexByGoTerm = new HashMap<String, Integer>(other.indexByGoTerm);
        this.geneById = new ArrayList<String>(other.geneById);
        this.goTermById = new ArrayList<String>(other.goTermById);
        this.values = new HashMap<Integer, Map<Integer, Double>>();

        for (int rowId : other.values.keySet()) {
            HashMap<Integer, Double> myRow = new HashMap<Integer, Double>(other.values.get(rowId));
            this.values.put(rowId, myRow);
        }
        this.countByGoTerm = new HashMap<Integer, Integer>(other.countByGoTerm);
    }

    public double getScoreForProteinAndGOterm(String protein, String GOTerm) {
        if (this.containsProtein(protein)) {
            int indexProt = this.indexByGene.get(protein);

            if (this.indexByGoTerm.containsKey(GOTerm)) {
                int indexGO = this.indexByGoTerm.get(GOTerm);

                if (this.values.get(indexProt).containsKey(indexGO)) {
                    return this.values.get(indexProt).get(indexGO);
                } else {
                    return 0.0;
                }

                //return this.values.get(indexProt).get(indexGO);
            } else {
                return 0.0;
            }
        } else {
            return 0.0;
        }

    }

    @Override
    public Assignment clone() {
        Assignment other = new Assignment(this);
        return other;
    }

    /**
     * adds a protein with an empty row
     */
    public void addProtein(final String proteinName) {
        if (!this.indexByGene.containsKey(proteinName)) {
            int id = this.geneById.size();
            this.geneById.add(proteinName);
            this.indexByGene.put(proteinName, id);
            this.values.put(id, new HashMap<Integer, Double>());
        }
    }

    public Assignment filterBySetOfGOTerms(final Set<String> validGOTerms) {
        HashSet<Integer> translatedValidIds = new HashSet<Integer>();

        HashSet<String> mineAndValid = new HashSet<String>(validGOTerms);
        mineAndValid.retainAll(this.indexByGoTerm.keySet());

        for (final String term : mineAndValid) {
            translatedValidIds.add(this.indexByGoTerm.get(term));
        }

        Assignment output = new Assignment();
        //int protId = 0;
        for (Map.Entry<Integer, Map<Integer, Double>> pairProtIdRow : this.values.entrySet()) {
            //output.addProtein(protein);
            String protein = this.geneById.get(pairProtIdRow.getKey());
            for (final int term_i : pairProtIdRow.getValue().keySet()) {

                if (translatedValidIds.contains(term_i)) {
                    final String goTermName = this.goTermById.get(term_i);
                    final double val = pairProtIdRow.getValue().get(term_i);
                    output.setValue(protein, goTermName, val);
                }
            }
        }

        return output;
    }

    public void setValue(final String protein, final String goTerm, final double value) {
        int idProt, idTerm;
        Map<Integer, Double> row;

        if (this.indexByGene.containsKey(protein)) {
            idProt = indexByGene.get(protein);
            row = this.values.get(idProt);
        } else {
            idProt = geneById.size();
            geneById.add(protein);
            indexByGene.put(protein, idProt);
            row = new HashMap<Integer, Double>();
            this.values.put(idProt, row);
        }

        if (this.indexByGoTerm.containsKey(goTerm)) {
            idTerm = indexByGoTerm.get(goTerm);
        } else {
            idTerm = goTermById.size();
            goTermById.add(goTerm);
            indexByGoTerm.put(goTerm, idTerm);
        }

        if (!row.containsKey(idTerm)) {
            int oldCount = countByGoTerm.containsKey(idTerm) ? countByGoTerm.get(idTerm) : 0;
            countByGoTerm.put(idTerm, oldCount + 1);
            row.put(idTerm, value);
        } else {
            // we only keep the maximum score, in case we try to overwrite the value
            row.put(idTerm, Math.max(value, row.get(idTerm)));
        }
    }

    public Assignment multiplyValuesByScalar(final double scalar) {

        Map<Integer, Map<Integer, Double>> values2 = new HashMap<Integer, Map<Integer, Double>>(values);

        for (int row : this.values.keySet()) {
            for (int col : this.values.get(row).keySet()) {
                values2.get(row).put(col, scalar * values.get(row).get(col));
            }
        }

        values = values2;
        return this;
    }

    public Set<String> getProteinsForGOTerm(final String goTermId) {
        Set<String> output = new HashSet<String>();

        if (this.indexByGoTerm.containsKey(goTermId)) {

            int id = this.indexByGoTerm.get(goTermId);

            for (Map.Entry<String, Integer> entry : indexByGene.entrySet()) {
                if (this.values.get(entry.getValue()).containsKey(id)) {
                    output.add(entry.getKey());
                }
            }
        }

        return output;
    }

    public Set<String> getRowIdentifiers() {
        return this.indexByGene.keySet();
    }

    public Map<Integer, Double> getRow(final String rowId) {
        return values.get(this.indexByGene.get(rowId));
    }

    public Set<String> getColumnIdentifiers() {
        return this.indexByGoTerm.keySet();
    }

    public Assignment add(Assignment other) {

        for (String row : other.getRowIdentifiers()) {

            int translatedIdRow = -1;
            Map<Integer, Double> targetRow = null;

            if (this.indexByGene.containsKey(row)) {
                translatedIdRow = this.indexByGene.get(row);
                targetRow = this.values.get(translatedIdRow);

            } else {
                translatedIdRow = this.geneById.size();
                this.geneById.add(row);
                this.indexByGene.put(row, translatedIdRow);
                targetRow = new HashMap<Integer, Double>();
                this.values.put(translatedIdRow, targetRow);
            }

            for (int untranslatedIdCol : other.getRow(row).keySet()) {
                double value = other.getRow(row).get(untranslatedIdCol);

                String otherColName = other.getGeneFromId(untranslatedIdCol);

                int translatedIdCol = -1;

                if (this.indexByGoTerm.containsKey(otherColName)) {
                    translatedIdCol = this.indexByGoTerm.get(otherColName);
                } else {
                    translatedIdCol = this.indexByGoTerm.size();
                    this.indexByGoTerm.put(otherColName, translatedIdCol);
                    this.goTermById.add(otherColName);
                }

                double newVal = value;

                if (targetRow.containsKey(translatedIdCol)) {
                    newVal += targetRow.get(translatedIdCol);
                }

                targetRow.put(translatedIdCol, newVal);
            }
        }

        return this;
    }

    /*
     public Matrix getMatrixAndGetRowColIds(List<String> protIds, List<String> goTermIds) {
     final int nrows = this.geneById.size();
     final int ncols = this.goTermById.size();
     Matrix matrix = new Matrix(nrows, ncols);

     for (int row : this.values.keySet()) {
     protIds.add(this.geneById.get(row));
     for (int col : this.values.get(row).keySet()) {
     matrix.set(row, col, this.values.get(row).get(col));

     if (goTermIds.size() != ncols) {
     goTermIds.add(this.goTermById.get(col));
     }
     }
     }

     return matrix;
     }

     public Matrix getMatrix(List<String> geneIds) {
     final int nrows = geneIds.size();
     final int ncols = this.goTermById.size();
     Matrix matrix = new Matrix(nrows, ncols);

     int rowCount = 0;

     for (String rowId : geneIds) {

     if (this.indexByGene.containsKey(rowId)) {
     int row = this.indexByGene.get(rowId);
     if (this.values.containsKey(row)) {
     for (int col : this.values.get(row).keySet()) {
     matrix.set(rowCount, col, this.values.get(row).get(col));
     }
     }
     }
     ++rowCount;
     }

     return matrix;
     }
    
     public Matrix getMatrix(Map<String, Integer> rowIds, Map<String, Integer> colIds)
     { 
     int M = rowIds.size();
     int N = colIds.size();
        
     Matrix ret = new Matrix(M,N);
        
     for (int ii : this.values.keySet())
     { 
     String s_ii = this.geneById.get(ii);
     int i = rowIds.get(s_ii);
            
     for (int jj : this.values.get(ii).keySet())
     {
     double val = this.values.get(ii).get(jj);
                
     String s_jj = this.goTermById.get(jj);
                
     int j = colIds.get(s_jj);
                
     ret.set(i, j, val);            
     }        
     }    
     return ret;
     }

     public Matrix getMatrix() {
     final int nrows = this.geneById.size();
     final int ncols = this.goTermById.size();
     Matrix matrix = new Matrix(nrows, ncols);

     for (int row : this.values.keySet()) {
     for (int col : this.values.get(row).keySet()) {
     matrix.set(row, col, this.values.get(row).get(col));
     }
     }

     return matrix;
     }
     */
    public int countNumberOfGenesForGOTerm(final String goTermId) {
        /*
         * if (this.indexByGoTerm.containsKey(goTermId)) { int termId =
         * this.indexByGoTerm.get(goTermId);
         *
         * for (Map<Integer, Double> mp : this.values.values()) { if
         * (mp.containsKey(termId)) { ++count; } } }
         */
        if (this.indexByGoTerm.containsKey(goTermId)) {
            int termId = this.indexByGoTerm.get(goTermId);
            return this.countByGoTerm.get(termId);
        } else {
            return 0;
        }

    }

    public Map<String, Double> getGOTermScoresForProteinId(final String proteinId) {
        Map<String, Double> result = new HashMap<String, Double>();

        for (Map.Entry<Integer, Double> pair : this.getRow(proteinId).entrySet()) {
            result.put(this.getGOTermFromId(pair.getKey()), pair.getValue());
        }
        return result;
    }

    public Assignment getFilteredVersionOnlyWithTermsAnnotatingBetweenMinMaxGenes(final int min, final int max) {

        // keep set of "allowed" terms
        Set<String> allowed = new HashSet<String>();
        for (String goTerm : this.goTermById) {
            int numGenes = this.countNumberOfGenesForGOTerm(goTerm);

            if (numGenes >= min && numGenes <= max) {
                allowed.add(goTerm);
            }
        }

        Assignment other = new Assignment();

        for (String proteinId : this.geneById) {
            Map<String, Double> row = this.getGOTermScoresForProteinId(proteinId);

            for (Map.Entry<String, Double> entry : row.entrySet()) {
                if (allowed.contains(entry.getKey())) {
                    other.setValue(proteinId, entry.getKey(), entry.getValue());
                }
            }
        }
        return other;
    }

    public int sizeGenes() {
        return this.geneById.size();
    }

    public int sizeTerms() {
        return this.goTermById.size();
    }

    public String getGeneFromId(final int untranslatedIdRow) {
        return this.geneById.get(untranslatedIdRow);
    }

    public String getGOTermFromId(final int untranslatedIdCol) {
        return this.goTermById.get(untranslatedIdCol);
    }

    public void readFromSeedFile(String inputFileName) throws IOException, NumberFormatException {
        BufferedReader reader = new BufferedReader(new FileReader(inputFileName));
        while (reader.ready()) {
            String line = reader.readLine();
            String fields[] = line.split("\\t");

            String protein = fields[0].trim();
            String goTerm = fields[1].trim();
            double score = Double.parseDouble(fields[2]);
            this.setValue(protein, goTerm, score);
        }
        reader.close();

    }

    public void writeSeedFile(String outputFileName) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputFileName));

        for (int protId : this.values.keySet()) {
            String protein = this.geneById.get(protId);
            for (int goTermId : this.values.get(protId).keySet()) {
                double score = this.values.get(protId).get(goTermId);
                String goTerm = this.goTermById.get(goTermId);

                writer.write(protein + "\t" + goTerm + "\t" + score);
                writer.newLine();
            }
        }

        writer.close();
    }
}
