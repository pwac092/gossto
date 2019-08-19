  ______   ______   ______   ______   ______  ______    
 /\  ___\ /\  __ \ /\  ___\ /\  ___\ /\__  _\/\  __ \   
 \ \ \__ \\ \ \/\ \\ \___  \\ \___  \\/_/\ \/\ \ \/\ \  
  \ \_____\\ \_____\\/\_____\\/\_____\  \ \_\ \ \_____\ 
    \/_____/ \/_____/ \/_____/ \/_____/   \/_/  \/_____/ 


# Introduction 

Welcome to GOssTo, the Gene Ontology Semantic Similarity Tool!

GOssTo is a software system for calculating semantic similarities between gene products in the Gene Ontology.
GOssTo implements the *Random Walk Contribution*, improving the accuracy of similarity measures.
For Yang's et. al description of the method:

* Improving GO semantic similarity measures by exploring the ontology beneath the terms and modelling uncertainty, *Bioinformatics*, **28**(10), pp. 1383-1389, 2012.

If you use GOssTo, please cite it.

* GOssTo: an extendable stand-alone and web tool for calculating semantic similarities on the Gene Ontology, *Bioinformatics*, **30**(15),  pp. 2235-2236, 2014.

Also, check the **Gossto paper website** (https://paccanarolab.org/gossto) to get a thorough documentation, and to download a bundle with the jar a several example dataset (Gene Ontology, etc.).

# Author labs' websites
These are the websites of the labs that contributed to Gossto.

* PaccanaroLab: http://www.paccanarolab.org
* Valentini Lab: https://sites.google.com/site/anacletolaboratory/

# Getting started

You can choose to run Gossto from the provided jar file (`dist/Gossto.jar`) or build it yourself. You can also download the jar file from the paper website.

## Building Gossto
###  Prerequisites

Gossto was designed without many dependencies. If you want to build it from sources you will need:

* Java (Java 5 or higher).
* Ant (version 1.9 or higher).

There is one additional dependency (Apache Commons Cli), which will be resolved by the use of Apache Ivy. In the Gossto source we also included a modified version of Jama (https://math.nist.gov/javanumerics/jama/), using single float precision, as this is allowed by its *public-domain license*.

### Using `ant` to build Gossto

In order to build Gossto, type 

```
ant jar
```

This will first fetch all dependencies, compile the sources, and then build `build/gossto.jar`. That jar file is a 'fat jar', meaning that the generated file will not have any external dependency, and then, it will be easy to run in any system where a Java VM is available.

## Running Gossto

Gossto can be run via the `Gossto.jar` file in two different manners:

* *Interactive way*: just run

```
java -jar Gossto.jar
```  
and follow the instructions.

* *Command-line way*: options can be added after the `Gossto.jar`, such as:
```
java -jar Gossto.jar -calculationdata termwise -calculationtype ism -evidencecodes EXP,IDA,IPI,IMP,IGI,IEP,TAS,IC -goapath gene association.goa yeast -obopath gene ontology ext.obo -hsm Resnik -hsmoutput demo hsm output -ismoutput demo ism output -ontology all -relations is a,part of -weightedJaccard true -terms all
```
A full list of options can be checked via the `java -jar Gossto.jar -help` option.

## Manual

There is a detailed manual on how to use Gossto, a FAQ, a list of possible errors, and several examples in the **Gossto paper website**: https://paccanarolab.org/gossto

## License

Gossto is *free software*, released under the GNU General Public License (GPL) v3.0. A copy of it is made available in the LICENSE file.

