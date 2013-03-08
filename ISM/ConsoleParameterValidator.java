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

import java.io.FileNotFoundException;
import java.io.IOException;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import util.TinyLogger;

/**
 *
 * @author aeromero
 */
public class ConsoleParameterValidator extends ParameterValidator {

    private final String[] args;

    public ConsoleParameterValidator(String args[]) {
        super();
        this.args = args;
    }

    @Override
    public void validate(IoValidation validate, TinyLogger logger) throws FileNotFoundException, IOException {
        CommandLineParser parser = new GnuParser();
        CommandLine cmd;
        HelpFormatter formatter = new HelpFormatter();
        //#####cli code#####
        //Setting up the parameter labels for console input
        Options paramOptions = new Options();
        
        paramOptions.addOption("help", false, "Shows this message");
        paramOptions.addOption("version", false, "v0.1 Aug 2012");
        paramOptions.addOption("getw", false, "Show warranty disclaimer section of the GNU GPL license");
        paramOptions.addOption("logfile", false, "Enter 'shallow' or 'deep' for the level of logging you require");
        paramOptions.addOption("obopath", true, "Enter filepath for OBO file");
        paramOptions.addOption("goapath", true, "Enter filepath for GOA file");
        paramOptions.addOption("relations", true, "Enter Gene Ontology relations to be used");
        paramOptions.addOption("evidencecodes", true, "Enter evidence codes to be used when parsing a GOA file");
        paramOptions.addOption("hsm", true, "Enter the name of the HSM to be used");
        paramOptions.addOption("ontology", true, "Enter 'bp' for Biological Process, 'mf' for Molecular Function, 'cc' for Cellular Component and 'all' for all three ontologies.");
        paramOptions.addOption("calculationtype", true, "Enter 'hsm' or 'ism' respective of the calculation you wish to undertake");
        paramOptions.addOption("calculationdata", true, "Enter 'termwise' or 'genewise' respective of the data you wish to base your calculation on");
        paramOptions.addOption("terms", true, "Enter the GO terms or gene IDs you wish to use");
        paramOptions.addOption("hsmoutput", true, "Enter the name and location the HSM outputs will be stored at");
        paramOptions.addOption("ismoutput", true, "Enter the name and location the ISM outputs will be stored at");
        paramOptions.addOption("weightedJaccard", true, "Select whether to weight the Jaccard Index for genewise ISM with the information content");
        
        try {
            cmd = parser.parse(paramOptions, args);
            //administrative parameters
            if (cmd.hasOption("help")) //prints list of parameter labels with description
            {
                formatter.printHelp("GOssTo", paramOptions);
                System.exit(-1);
            }
            if (cmd.hasOption("version")) //prints GOssTo version
            {
                System.out.println("Version : " + paramOptions.getOption("version").getDescription());
                System.exit(-1);
            }
            if (cmd.hasOption("getw")) //Prints the GNU GPL warranty section
            {
                System.out.println("15. Disclaimer of Warranty. \n THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. \n EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR \n OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY KIND,  \n EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES \n OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE ENTIRE RISK \n AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM \n PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.");
                System.exit(-1);
            }
            //this.log file parameter
            if (cmd.hasOption("this.logfile")) //initialises this.log file
            {
                TinyLogger.setLogging(true);
                //this.log = true;
                logger.initialiseLogWriter();
            }
            //data parameters
            if (cmd.hasOption("obopath")) //OBO file path
            {
                this.oboFile = cmd.getOptionValue("obopath");
                this.notes.add(this.oboFile);
            } else {
                logger.logAndCloseWriter("############ ERROR: OBOpath not spcified");
                System.err.println("ERROR: OBO Path Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("goapath")) //GOA file path
            {
                this.goaFile = cmd.getOptionValue("goapath");
                this.notes.add(this.goaFile);
            } else {
                logger.logAndCloseWriter("############ ERROR: GOA path not specified");
                System.err.println("ERROR: GOA Path Not Specified");
                System.exit(-1);
            }
            validate.validateFilePaths(this.oboFile, this.goaFile);
            if (cmd.hasOption("relations")) //GO Relations
            {
                String rawRelations = cmd.getOptionValue("relations");
                this.chosenRelations = rawRelations.split(",");
                String temp = "";
                for (int i = 0; i < this.chosenRelations.length; i++) {
                    this.chosenRelations[i] = chosenRelations[i].toLowerCase();
                    temp += this.chosenRelations[i] + " ";
                }
                this.notes.add(temp);
                validate.validateRelations(this.chosenRelations);
            } else {
                logger.logAndCloseWriter("############ ERROR: Go relations not specified");
                System.err.println("ERROR: GO Relations Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("evidencecodes")) //Evidence Codes
            {
                String rawCodes = cmd.getOptionValue("evidencecodes");
                this.evidenceCodes = rawCodes.split(",");
                String temp = "";
                for (int i = 0; i < this.evidenceCodes.length; i++) {
                    this.evidenceCodes[i] = this.evidenceCodes[i].toUpperCase(); //needed as upper case for parsing
                    temp += this.evidenceCodes[i] + " ";
                }
                this.notes.add(temp);
                validate.validateEvidenceCodes(this.evidenceCodes);
            } else {
                logger.logAndCloseWriter("############ ERROR: Evidence codes not specified");
                System.err.println("ERROR: Evidence Codes Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("hsm")) //HSM method
            {
                this.hsmChoice = cmd.getOptionValue("hsm");
                this.notes.add(this.hsmChoice);
            } else {
                logger.logAndCloseWriter("############ ERROR: HSM not specified");
                System.err.println("ERROR: HSM Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("ontology")) //GO Ontology choice
            {
                this.dagChoice = cmd.getOptionValue("ontology").toLowerCase();
                this.dagChoice = validate.validateDagChoice(dagChoice);
                this.notes.add(this.dagChoice);
            } else {
                logger.logAndCloseWriter("############ ERROR: Ontology not specified");
                System.err.println("ERROR: Ontology Choice Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("calculationtype")) //Whether calculating HSM or ISM
            {
                if (cmd.getOptionValue("calculationtype").toLowerCase().equals("ism") == true) {
                    this.ismChoice = true;
                    this.notes.add("ISM");
                } else if (cmd.getOptionValue("calculationtype").toLowerCase().equals("hsm") == true) {
                    this.ismChoice = false;
                    this.notes.add("HSM");
                } else {
                    logger.logAndCloseWriter("############ ERROR: Calculation type entered improperly");
                    System.err.println("ERROR: HSM/ISM specification has been entered incorrectly");
                }
            } else {
                logger.logAndCloseWriter("############ ERROR: Calculation type not set to HSM  or ISM");
                System.err.println("ERROR: It Has Not Been Specified Whether an HSM or ISM is Being Calculated");
                System.exit(-1);
            }
            if (cmd.hasOption("calculationdata")) //Whether we're working with gene or GO term data
            {
                if (cmd.getOptionValue("calculationdata").toLowerCase().equals("genewise")) {
                    this.termWise = false;
                    this.notes.add("genewise");
                } else if (cmd.getOptionValue("calculationdata").toLowerCase().equals("termwise")) {
                    this.termWise = true;
                    this.notes.add("termwise");
                } else {
                    logger.logAndCloseWriter("############ ERROR: termwise or genewise choice entered improperly");
                    System.err.println("ERROR: Genewise/Termwise specification has been entered incorrectly");
                    System.exit(-1);
                }
            } else {
                logger.logAndCloseWriter("############ ERROR: termwise or genewise not set");
                System.err.println("ERROR: It Has Not Been Specified Whether Gene or GO Terms are to be Used");
                System.exit(-1);
            }
            if (cmd.hasOption("terms")) //Specific GO terms or genes to be used
            {
                if (!cmd.getOptionValue("terms").toLowerCase().equals("all")) {
                    String rawTerms = cmd.getOptionValue("terms");
                    if (this.termWise == true) {
                        this.goIDs = rawTerms.split(",");
                    } else {
                        this.geneIDs = rawTerms.split(",");
                    }
                    this.notes.add("Terms used: " + cmd.getOptionValue("terms"));
                } else {
                    if (this.termWise == true) {
                        this.notes.add("all GO terms");
                    } else {
                        this.notes.add("all genes");
                    }
                }
            } else {
                logger.logAndCloseWriter("############ ERROR: Terms parameter not set properly");
                System.err.println("ERROR: Terms (Gene IDs or GO Terms) Not Specified or Not Set To 'all'");
                System.exit(-1);
            }
            if (cmd.hasOption("hsmoutput")) //HSM output file path & name
            {
                this.hsmFileName = cmd.getOptionValue("hsmoutput");
                validate.validateOutputLocation(this.hsmFileName);
            } else {
                logger.logAndCloseWriter("############ ERROR: No output path for HSM");
                System.err.println("ERROR: HSM Output Path Not Specified");
                System.exit(-1);
            }


            if (cmd.hasOption("ismoutput")) //ISM output file path and name
            {
                this.ismFileName = cmd.getOptionValue("ismoutput");
                validate.validateOutputLocation(this.ismFileName);
            } else if (this.ismChoice == true) {
                logger.logAndCloseWriter("############ ERROR: No output path for ISM");
                System.err.println("ERROR: ISM Output Path Not Specified");
                System.exit(-1);
            }
            if (cmd.hasOption("weightedJaccard")) //check for weightedJaccard
            {
                if (cmd.getOptionValue("weightedJaccard").toLowerCase().equals("true")) {
                    this.weightedJaccard = true;
                } else if (cmd.getOptionValue("weightedJaccard").toLowerCase().equals("false")) {
                    this.weightedJaccard = false;
                } else { //just in case something weird was written.
                    logger.logAndCloseWriter("############ ERROR: Invalid choice for weightedJaccard option");
                    System.err.println("ERROR: Invalid choice for weightedJaccard option");
                    System.exit(-1);
                }
            } else if (this.ismChoice == true) { //quick check to verify if the option is set 
                logger.logAndCloseWriter("############ ERROR: weightedJaccard option has to be set when computing ISM");
                System.err.println("ERROR: weightedJaccard option has to be set when computing ISM");
                System.exit(-1);
            }

        } catch (ParseException e) {
            logger.logAndCloseWriter("############ ERROR: Parse Failed");
            System.err.println("ERROR: Parse failed : " + e.getMessage());
            System.exit(-1);
        }
    }
}
