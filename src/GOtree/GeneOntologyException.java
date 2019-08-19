/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package GOtree;

/**
 *
 * @author aeromero
 */
public class GeneOntologyException extends Exception {
    String my_message;

    public GeneOntologyException(String message) {
        super(message);
        this.my_message = message;
    }

    public String getMy_message() {
        return my_message;
    }    
}
