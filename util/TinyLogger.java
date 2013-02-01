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
package util;

import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

/*
 * A small and not-so-powerful logger
 * @author Alfonso E. Romero
 */
public class TinyLogger {

    private BufferedWriter logwriter; //Used for writing messages to the log file
    private static boolean logging = false;

    public TinyLogger() {
    }

    public static void setLogging(boolean _logging) {
        logging = _logging;
    }

    public void initialiseLogWriter(String fileName) throws FileNotFoundException, IOException {
        logwriter = new BufferedWriter(new FileWriter(fileName));
    }

    //Initialises the writer for the log file
    public void initialiseLogWriter() throws FileNotFoundException, IOException {
        this.initialiseLogWriter("LOG.txt");
    }

    public void log(String log) throws IOException {
        if (logging) {
            logwriter.write(log);
            logwriter.newLine();
        }
    }

    public void logAndCloseWriter(String log) throws IOException {
        if (logging) {
            log(log);
            this.logwriter.close();
        }
    }

    public void showMessage(String message) throws IOException {
        System.out.println(message);
        this.log(message);
    }
}
