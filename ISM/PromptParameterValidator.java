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

import java.io.FileNotFoundException;
import java.io.IOException;
import util.TinyLogger;

/**
 *
 * @author aeromero
 */
public class PromptParameterValidator extends ParameterValidator {

    private TerminalInterface ui; //interface object interacts with the user

    public PromptParameterValidator(TerminalInterface ui_) {
        super();
        ui = ui_;
    }

    @Override
    public void validate(IoValidation validate, TinyLogger logger) throws FileNotFoundException, IOException {
        this.oboFile = ui.requestOBOfile(); //OBO file path
        this.notes.add(this.oboFile);
        this.goaFile = ui.requestGOAfile(); //GOA file path
        this.notes.add(this.goaFile);
        validate.validateFilePaths(this.oboFile, this.goaFile);
        this.chosenRelations = ui.requestRelations(); //GO Relations
        String temp = "";
        for (String rel : this.chosenRelations) {
            temp += rel + " ";
        }
        this.notes.add(temp);
        validate.validateRelations(this.chosenRelations);
        this.evidenceCodes = ui.requestEvidenceCodes(); //Evidence Codes
        for (int i = 0; i < this.evidenceCodes.length; i++) {
            this.evidenceCodes[i] = this.evidenceCodes[i].toUpperCase(); //needed as upper case for parsing
        }
        validate.validateEvidenceCodes(this.evidenceCodes);
        temp = "";
        for (String eviC : this.evidenceCodes) {
            temp += eviC + " ";
        }
        this.notes.add(temp);

        this.termWise = ui.requestTermwise(); //Whether Genewise or Termwise
        if (ui.usingOwnHSM() == true) //Whether to pick HSM from list or to type in the name
        {
            this.hsmChoice = ui.requestOwnHSMname();
            if (this.hsmChoice.equals("Resnik") == true || hsmChoice.equals("Jiang") == true || hsmChoice.equals("Lin") == true || hsmChoice.equals("simGraSM") == true || hsmChoice.equals("simUI") == true || hsmChoice.equals("simGIC") == true) {
                this.ownHSM = false;
            } else {
                this.ownHSM = true; //if name is not one of the 6 supplied HSMs
            }
            if (this.termWise == false && this.ownHSM == true) {
                this.weightedJaccard = ui.requestJaccardChoice(); //Ask user whether they'd like to use a weighted jaccard index
            }
        } else {
            this.hsmChoice = ui.requestHSMChoice(this.termWise); //list the 6 supplied HSM choices for selection
        }
        this.notes.add(this.hsmChoice);

        //GET THE DESIRED ONTOLOGY
        this.dagChoice = validate.validateDagChoice(ui.requestDAG()); //Choice of GO Ontology
        this.notes.add(this.dagChoice);
        if (this.termWise == false) {
            this.geneIDs = ui.requestGeneIDs(); //Request any specific Gene IDs
        } else {
            this.goIDs = ui.requestCategories(); //Request any specific GO IDs
        }
        //GET WHETHER WE CALCULATE ISM OR NOT.
        this.ismChoice = ui.requestISMChoice(); //Request whether we are calculating an HSM or an ISM
        if (this.ismChoice == false) {
            this.notes.add("HSM");
            this.hsmFileName = ui.requestFileName(true); //Request HSM name & file path
            this.ismFileName = "";
        } else {
            this.notes.add("ISM");
            this.hsmFileName = ui.requestFileName(true); //Request HSM name & file path
            validate.validateOutputLocation(this.hsmFileName);
            this.ismFileName = ui.requestFileName(false); //Request ISM name & file path
            validate.validateOutputLocation(this.ismFileName);
        }

        if (this.termWise == true) {
            this.notes.add("termwise");
            if (this.goIDs != null) {
                String gotemp = "Terms used: ";
                for (String go : this.goIDs) {
                    gotemp += go + ",";
                }
                this.notes.add(gotemp);
            } else {
                this.notes.add("all GO terms");
            }
        } else {
            this.notes.add("genewise");
            if (this.geneIDs != null) {
                String genetemp = "Terms used: ";
                for (String gene : this.geneIDs) {
                    genetemp += gene + ",";
                }
                this.notes.add(genetemp);
            } else {
                this.notes.add("all genes");
            }
        }
    }
}
