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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 *
 * @author Samuel Heron
 */
//Deals with user interaction. Prompts user for input parameters
public class TerminalInterface {

    private Scanner userInput; //reads user inputed values from the command line / console

    //instantiates the scanner variable
    TerminalInterface() {
        userInput = new Scanner(System.in); //initialise scanner to read from the console
    }

    //Requests a file path for the Gene Ontology OBO file
    public String requestOBOfile() {
        System.out.println("Please enter the file path for the Gene Ontology (OBO) file you wish to use: ");
        System.out.print(">? ");
        String OBOpath;
        OBOpath = userInput.next();
        return OBOpath;
    }

    //Requests a file path for the GOA Annotation file
    public String requestGOAfile() {
        System.out.println("Please enter the file path for the Annotation (GOA) file you wish to use: ");
        System.out.print(">? ");
        String GOApath;
        GOApath = userInput.next();
        return GOApath;
    }

    //Requests a file name & location for the similarity results, the parameter decides which file/file group
    public String requestFileName(boolean HSM) {
        String name = "", measure = "";
        if (HSM == true) {
            measure = "HSM";
        } else {
            measure = "ISM";
        }
        System.out.println("Please enter the desired file name for the " + measure + " results:");
        System.out.print(">? ");
        name = userInput.next();
        return name;
    }

    //Asks the user whether they'd like to pick the HSM method from a list or enter its name
    public boolean usingOwnHSM() {
        String choice = "";
        System.out.println("Would you like to pick an HSM from a list? (y) or enter the name of the desired HSM? (n)");
        System.out.print(">? ");
        choice = this.userInput.next();
        if (choice.toLowerCase().equals("y") == true) {
            return false;
        } else if (choice.toLowerCase().equals("n") == true) {
            return true;
        } else {
            incorrect(choice);
            return false; //to make the debugger happy
        }
    }

    //Error message for incorrect parameter entry/choice by the user
    private void incorrect(String choice) {
        System.err.println("ERROR: " + choice + " is not a valid option");
        System.exit(0);
    }

    //request whether the user would like to use a weighted jaccard index
    public boolean requestJaccardChoice() {
        String choice = "";
        System.out.println("Does your HSM use a weighted Jaccard Index in its calculations?: (y/n)");
        System.out.print(">? ");
        choice = this.userInput.next();
        if (choice.toLowerCase().equals("y") == true) {
            return true;
        } else if (choice.toLowerCase().equals("n") == true) {
            return false;
        } else {
            incorrect(choice);
            return false; //to make the debugger happy
        }
    }

    //Requests the HSM to be calculated and used by the ISM, the parameter ensures only relevant HSMs are listed
    public String requestHSMChoice(boolean termwise) {
        String chosenHSM = "";
        int hsmChoice = 0;
        System.out.println("Please select the Host Similarity Measure (HSM) that you would like to use from the following list: ");
        System.out.println();
        System.out.println("Termwise & Genewise Measures:");
        System.out.println("1. Resnik");
        System.out.println("2. Lin");
        System.out.println("3. Jiang");
        System.out.println("4. simGraSM");
        System.out.println("");
        if (termwise == false) {
            System.out.println("Solely Genewise Measures:");
            System.out.println("5. simUI");
            System.out.println("6. simGIC");
            System.out.println("");
        }
        System.out.println("Please make your choice by entering the number of the desired measure.");
        System.out.print(">? ");
        hsmChoice = this.userInput.nextInt();
        switch (hsmChoice) {
            case 1: {
                chosenHSM = "Resnik";
                break;
            }
            case 2: {
                chosenHSM = "Lin";
                break;
            }
            case 3: {
                chosenHSM = "Jiang";
                break;
            }
            case 4: {
                chosenHSM = "simGraSM";
                break;
            }
            case 5: {
                if (termwise == true) {
                    incorrect(Integer.toString(hsmChoice));
                }
                chosenHSM = "simUI";
                break;
            }
            case 6: {
                if (termwise == true) {
                    incorrect(Integer.toString(hsmChoice));
                }
                chosenHSM = "simGIC";
                break;
            }
            default: {
                incorrect(Integer.toString(hsmChoice));
                System.exit(0);
                break;
            }
        }
        return chosenHSM;
    }

    //Get name of user supplied HSM
    public String requestOwnHSMname() {
        System.out.println("Please enter the name of the HSM you wish to use:");
        System.out.print(">? ");
        String desired_class_name = this.userInput.next();
        return desired_class_name;
    }

    //Request whether to use Genes or GO terms as a basis for the calculation
    public boolean requestTermwise() {
        String choice = "";
        System.out.println("Would you like to get genewise (g) or termwise (t) results?: (g / t)");
        System.out.print(">? ");
        choice = this.userInput.next();
        if (choice.toLowerCase().equals("g") == true) {
            return false;
        } else if (choice.toLowerCase().equals("t") == true) {
            return true;
        } else {
            incorrect(choice);
            return false; //to make the debugger happy
        }
    }

    //Asks the user whether or not they would like to calculate the ISM or just the HSMs
    public boolean requestISMChoice() {
        boolean ismChoice = true;
        String choice = "";
        System.out.println("Would you like to calculate solely HSM results? (y / n)");
        System.out.print(">? ");
        choice = this.userInput.next();
        if (choice.toLowerCase().equals("y") == true) {
            ismChoice = false;
        } else if (choice.toLowerCase().equals("n") == true) {
            ismChoice = true;
        } else {
            incorrect(choice);
        }
        return ismChoice;
    }

    public String[] requestEvidenceCodes() {
        String[] choice = null;
        System.out.println("Would you like to choose which evidence codes will be accepted? (y / n for default)");
        System.out.println("The default list is: EXP, IDA, IPI, IMP, IGI, IEP, TAS, IC"); //default
        System.out.print(">? ");
        String yesOrNo = this.userInput.next();

        if (yesOrNo.toLowerCase().equals("y") == true) {
            boolean isAProblem;
            int count = 0;
            String[] available = new String[]{"exp", "ida", "ipi", "imp", "igi", "iep", "iss", "iso", "isa", "ism", "igc", "iba", "ibd", "ikr", "ird", "rca", "tas", "nas", "ic", "nd", "iea"};

            HashMap<String, Integer> countAppearances = new HashMap<String, Integer>();

            do {
                isAProblem = false;

                System.out.println("The available evidence codes are:");
                //just to print the available evidence codes nicely.
                for (int i = 0; i < available.length; i++) {
                    System.out.print(available[i].toUpperCase());
                    if (available.length - 1 > i) {
                        System.out.print(",");
                    }
                }
                System.out.println();
                //---

                for (String rel : available) {
                    countAppearances.put(rel, 0);
                }

                System.out.println("List the evidence codes separating them with a comma. Press enter once you are done.");
                System.out.print(">? ");

                String providedList = this.userInput.next().toLowerCase();
                for (String current : providedList.split(",")) {
                    current = current.trim();
                    if (!countAppearances.containsKey(current)) {
                        System.out.println("The specified evidence code " + current + " is not valid");
                        isAProblem = true;
                        break;
                    } else if (countAppearances.get(current) > 0) {
                        System.out.println("The specified evidence code " + current + " is repeated several times");
                        isAProblem = true;
                        break;
                    }
                    countAppearances.put(current, 1);
                    count++;
                }

            } while (isAProblem);
            String evCodes[] = new String[count];
            int index = 0;
            for (String ev : countAppearances.keySet()) {
                if (countAppearances.get(ev) > 0) {
                    evCodes[index] = ev;
                    ++index;
                }
            }
            return evCodes;
        } else //use the default.
        {
            choice = new String[]{"EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "TAS", "IC"};
        }
        return choice;
    }

    //Asks the user which ontology they would like to use
    public String requestDAG() {
        int choice = 0;

        System.out.println("Please select your choice of ontology from the list below by entering its corresponding number: ");
        System.out.println();
        System.out.println("1. Biological Process");
        System.out.println("2. Molecular Function");
        System.out.println("3. Cellular Component");
        System.out.println("4. All Three");
        System.out.print(">? ");
        choice = userInput.nextInt();

        switch (choice) {
            case 1:
                return "bp";
            case 2:
                return "mf";
            case 3:
                return "cc";
            case 4:
                return "all";
            default:
                incorrect(Integer.toString(choice));
                return "incorrect"; //to make the debugger happy
        }
    }

    //Requests any GO ID's to be worked with
    public String[] requestCategories() {


        System.out.println("Would you like to work with a specific set of GO terms? (y/n)");
        System.out.println("By selecting 'n' all GO terms will be considered");
        System.out.print(">? ");
        String choice = userInput.next();

        if (choice.toLowerCase().equals("y")) {
            //we will take in the goids input by the user. If a goterm is duplicated 
            //we will just ignore it. 
            //the format will be: go:0000009
            //we'll check the validity of the go term with a regular expression. nothing fancy
            Set<String> selectedGOCategories = new HashSet<String>();
            boolean success = true;
            do {
                System.out.println("You will be able to add as many GOterms as you wish. "
                        + ".\n Input the GOterms as a comma-separated value list.\n "
                        + "The format for the GOterms: GO:0000000 \n"
                        + "Repeated GOterms will be ignored. \n");
                //this will store the goterms and make sure there are no repeats
                //this will read the single line with all the goterms
                String providedList = this.userInput.next().toLowerCase();
                //now go through them one at a time.
                //so..match GO: first. after that, convert the rest of the string 
                //into an atomic unit, and check for any digit. find exactly seven of them.
                //String goTermPattern = "go:\\d.*/{7,/}";
                String goTermPattern = "go:\\d{7}";
                success = true;
                for (String current : providedList.split(",")) {
                    //1. check the format of "current"
                    if (!current.trim().toLowerCase().matches(goTermPattern)) {
                        System.out.println("The provided category string " + current + " does not match a valid GO category. Check the format.");
                        success = false;
                        break;
                    }
                    selectedGOCategories.add(current.trim().toLowerCase());
                }
            } while (!success);
            //just in case the user wants to troll us. 
            if (selectedGOCategories.isEmpty()) {
                return null;
            }
            //return the actual list.
            return selectedGOCategories.toArray(new String[selectedGOCategories.size()]);
        } else if (!choice.toLowerCase().equals("n")) {
            incorrect(choice);
        }
        //this function returns null to indicate thall all goterms are going to be used.
        return null;
    }

//Requests the IDs of particular genes the user would like to calculate the similarity of
    public String[] requestGeneIDs() {
        String[] GeneIDs = null;
        String choice = "";
        System.out.println("Would you like to work with a specific set of genes? (y / n)");
        System.out.print(">? ");
        choice = userInput.next();
        if (choice.toLowerCase().equals("y")) {
            int noGenes = 0;
            System.out.println("Please enter the number of Genes you wish to work with: ");
            System.out.print(">? ");
            noGenes = userInput.nextInt();
            GeneIDs = new String[noGenes];
            System.out.println("Please enter the Gene ID's you wish to work with one by one: ");
            for (int i = 0; i < noGenes; i++) {
                System.out.print(">? ");
                GeneIDs[i] = userInput.next();
            }
            return GeneIDs;
        } else if (choice.toLowerCase().equals("n")) {
            return GeneIDs;
        } else {
            incorrect(choice);
            return GeneIDs; //keep the debugger happy
        }
    }

    //Requests the relations the user wishes to use
    public String[] requestRelations() {
        String choice = "";
        System.out.println("By default this program will only use the 'is_a' GO relationship, would you like to include more relations? (y / n for default)");
        System.out.print(">? ");
        choice = userInput.next();

        if (choice.toLowerCase().equals("y") == true) {
            boolean isAProblem;
            int count = 0;
            String availableRelations[] = {"is_a", "part_of", "regulates", "positively_regulates", "negatively_regulates", "has_part"};

            HashMap<String, Integer> countAppearances = new HashMap<String, Integer>();

            do {
                isAProblem = false;

                System.out.println("The available relations codes are:");
                //just to print the available evidence codes nicely.
                for (int i = 0; i < availableRelations.length; i++) {
                    System.out.print(availableRelations[i].toUpperCase());
                    if (availableRelations.length - 1 > i) {
                        System.out.print(",");
                    }
                }
                System.out.println();
                //---


                for (String rel : availableRelations) {
                    countAppearances.put(rel, 0);
                }

                System.out.println("List the relations separating them with a comma. Press enter once you are done.");
                System.out.println("Keep in mind that the 'is_a' relation is always added.");
                System.out.print(">? ");

                String providedRelations = this.userInput.next().toLowerCase();
                for (String providedRelation : providedRelations.split(",")) {
                    providedRelation = providedRelation.trim();
                    if (!countAppearances.containsKey(providedRelation)) {
                        System.out.println("The specified relation " + providedRelation + " is not valid");
                        isAProblem = true;
                        break;
                    } else if (countAppearances.get(providedRelation) > 0) {
                        System.out.println("The specified relation " + providedRelation + " is repeated several times");
                        isAProblem = true;
                        break;
                    }
                    countAppearances.put(providedRelation, 1);
                    count++;
                }

            } while (isAProblem);

            String relations[] = new String[count];
            int index = 0;
            for (String rel : countAppearances.keySet()) {
                if (countAppearances.get(rel) > 0) {
                    relations[index] = rel;
                    ++index;
                }
            }
            return relations;
        } else if (choice.toLowerCase().equals("n") == true) {
            String relations[] = new String[]{"is_a"}; //default
            return relations;
        } else {
            incorrect(choice);
            return null; //keep the debugger happy 
        }
    }

    public void farewell() {

        System.out.println("#####Printing Complete#####");
        System.out.println();
        System.out.println("Thank you for using GOssTo.");
    }

    //welcome message displayed on launch of GOssTo, contains GNU GPL info
    public void welcomer() {
        System.out.println();
        System.out.println("******************************************************************************************************************");

        System.out.println("Welcome to GOssTo the Gene Ontology Semantic Similarity Tool!");
        System.out.print("GOssTo is a software system for ");
        System.out.println("calculating semantic similarities between gene products in the Gene Ontology.");
        System.out.print("GOssTo implements the Random Walk Contribution, improving the accuracy ");
        System.out.print("of similarity measures.");

        System.out.println();
        System.out.println("For Yang's et. al description of the method:");
        System.out.println();
        System.out.println("Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty");
        System.out.println("Bioinformatics, vol. 28, iss. 10, pp. 1383-1389, 2012.");

        System.out.println();
        System.out.println("If you use GOssTo, please cite it.");
        System.out.println();
        System.out.println("GOssTo: an extendable stand-alone and web tool for calculating semantic similarities on the Gene Ontology");
        System.out.println("To appear");

        System.out.println();
        System.out.println("PaccanaroLab: http://www.paccanarolab.org");
        System.out.println("Valentini Lab: http://homes.di.unimi.it/~valenti/index.html");
        System.out.println("********************************************************************************************************************");
        System.out.println("--------------------------------------------------------------------------------------------------------------------");
        System.out.println("GOssTo  Copyright (C) 2012  Samuel Heron, Alfonso E. Romero");
        System.out.println("This program comes with ABSOLUTELY NO WARRANTY; for details type `--getw' after the Jar execution statement.");
        System.out.println("This is free software, and you are welcome to redistribute it");
        System.out.println("under certain conditions; see the full GNU GPL licence for details.");
        System.out.println("The license can be found in full here: http://www.gnu.org/licenses/gpl.html");
        System.out.println("--------------------------------------------------------------------------------------------------------------------");
        System.out.println();
    }
}
