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

import java.util.Set;

/**
 *
 * @author aeromero
 */
public class AnnotationFileWithoutExternalProteinIds implements AnnotationFileStrategy {

    public AnnotationFileWithoutExternalProteinIds() {
    }

    @Override
    public void setProteinIdentifiers(Set<String> ids) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Set<String> getSetOfUnParsedProteins() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean proteinIdIsValid(String[] fieldsGOAnnotation) {
        return true;
    }

    @Override
    public String getProteinId(String[] fieldsGOAnnotation) {
        return fieldsGOAnnotation[2];
    }
}
