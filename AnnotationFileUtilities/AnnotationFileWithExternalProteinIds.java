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

/*
 * To change this template, choose Tools | Templates
 * and open the template in th e editor.
 */
package AnnotationFileUtilities;

import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author aeromero
 */
public class AnnotationFileWithExternalProteinIds implements AnnotationFileStrategy {

    Set<String> proteinIdentifiers;
    Set<String> parsedProteinIdentifiers;
    String validIdentifier;

    public AnnotationFileWithExternalProteinIds() {
        this.parsedProteinIdentifiers = new HashSet<String>();
    }

    @Override
    public void setProteinIdentifiers(Set<String> ids) {
        this.proteinIdentifiers = ids;
    }

    @Override
    public Set<String> getSetOfUnParsedProteins() {

        HashSet<String> unparsedProteinIdentifiers = new HashSet<String>(this.proteinIdentifiers);
        unparsedProteinIdentifiers.removeAll(this.parsedProteinIdentifiers);
        return unparsedProteinIdentifiers;
    }

    /*
    public String getProteinByField(String[] fields) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
*/
    @Override
    public boolean proteinIdIsValid(String[] fieldsGOAnnotation) {

        if (this.proteinIdentifiers.contains(fieldsGOAnnotation[1])) {
            this.validIdentifier = fieldsGOAnnotation[1];
            this.parsedProteinIdentifiers.add(validIdentifier);
            return true;
        } else if (this.proteinIdentifiers.contains(fieldsGOAnnotation[9])) {
            validIdentifier = fieldsGOAnnotation[9];
            this.parsedProteinIdentifiers.add(validIdentifier);
            return true;
        } else {
            String aliases[] = fieldsGOAnnotation[10].split("|");
            for (String alias : aliases) {
                if (this.proteinIdentifiers.contains(alias)) {
                    this.validIdentifier = alias;
                    this.parsedProteinIdentifiers.add(validIdentifier);
                    return true;
                }
            }
        }

        return false;
    }

    @Override
    public String getProteinId(String[] fieldsGOAnnotation) {
        return this.validIdentifier;
    }
}
