// ******************************************************************************
//
// Title:       Force Field X.
// Description: Force Field X - Software for Molecular Biophysics.
// Copyright:   Copyright (c) Michael J. Schnieders 2001-2025.
//
// This file is part of Force Field X.
//
// Force Field X is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License version 3 as published by
// the Free Software Foundation.
//
// Force Field X is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// Force Field X; if not, write to the Free Software Foundation, Inc., 59 Temple
// Place, Suite 330, Boston, MA 02111-1307 USA
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the
// GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a
// module which is not derived from or based on this library. If you modify this
// library, you may extend this exception to your version of the library, but
// you are not obligated to do so. If you do not wish to do so, delete this
// exception statement from your version.
//
// ******************************************************************************
package ffx.numerics.integrate;

import java.text.DecimalFormat;

import static org.apache.commons.math3.util.FastMath.abs;
import static org.apache.commons.math3.util.FastMath.cos;
import static org.apache.commons.math3.util.FastMath.sin;

/**
 * This program integrates using three methods: the trapezoidal method, Simpson's Three Point
 * Integration, and Boole's Five Point Integration
 *
 * @author Claire O'Connell
 */
public class Integration {

  private static final double[] x = new double[202];

  static {
    x[0] = 0;
    for (int i = 1; i < 201; i++) {
      x[i] = .0025 + .005 * (i - 1);
    }
    x[201] = 1;
  }

  private Integration() {
    // Prevent instantiation.
  }


  /**
   * averageIntegral.
   *
   * @param leftInt  a double.
   * @param rightInt a double.
   * @return a double.
   */
  public static double averageIntegral(double leftInt, double rightInt) {
    return (leftInt + rightInt) / 2.0;
  }

  /**
   * generateTestData_v1.
   *
   * @return an array of double values.
   */
  public static double[] generateTestData_v1() {
    double[] y = new double[202];

    for (int i = 0; i < 202; i++) {
      y[i] = 10 * sin(6 * x[i]) - 7 * cos(5 * x[i]) + 11 * sin(8 * x[i]);
    }

    return y;
  }

  /**
   * halfBinComposite.
   *
   * @param inputData an array of double values.
   * @param mode      the integration mode.
   * @param side      a {@link java.lang.String} object.
   * @return a double.
   */
  public static double halfBinComposite(double[] inputData, int mode, String side) {
    int n = inputData.length;
    double halfBinComposite = 0;
    // Split by side first, then integration type
    if (side.equalsIgnoreCase("left")) {
      // Using trapezoidal integral for lower half bin
      double lowerHalfBin = (inputData[1] + inputData[0]) / 2.0 * (x[1] - x[0]);
      halfBinComposite += lowerHalfBin;
      switch (mode) {
        // Case 1 is the Simpson's method, uses trapezoid on right for bin left out of Simpson's
        case 1 -> {
          double upperTrapArea = (inputData[n - 3] + inputData[n - 2]) / 2.0 * (x[n - 2] - x[n - 3]);
          halfBinComposite += upperTrapArea;
        }
        // Case 2 is the Boole's method, uses Simpsons and trapezoidal integral on right to cover
        // remaining bins
        case 2 -> {
          double upperSimpson =
              (1.0 / 3.0)
                  * (x[n - 4] - x[n - 5])
                  * (inputData[n - 5] + 4 * inputData[n - 4] + inputData[n - 3]);
          halfBinComposite += upperSimpson;
          double upperTrapArea = (inputData[n - 3] + inputData[n - 2]) / 2.0 * (x[n - 2] - x[n - 3]);
          halfBinComposite += upperTrapArea;
        }
      }
    } else if (side.equalsIgnoreCase("right")) {
      // Upper half bin calculated with trapezoid
      double upperHalfBin = (inputData[n - 1] + inputData[n - 2]) / 2.0 * (x[n - 1] - x[n - 2]);
      halfBinComposite += upperHalfBin;
      switch (mode) {
        // Case 1 is the Simpson's method, uses trapezoid on left for bin left out of Simpson's
        case 1:
          double lowerTrapArea = (inputData[1] + inputData[2]) / 2.0 * (x[2] - x[1]);
          halfBinComposite += lowerTrapArea;
          break;
        // Case 2 is the Boole's method, uses Simpsons and trapezoidal integral on left to cover
        // remaining bins
        case 2:
          lowerTrapArea = (inputData[1] + inputData[2]) / 2.0 * (x[2] - x[1]);
          halfBinComposite += lowerTrapArea;
          double lowerSimpson =
              (1.0 / 3.0) * (x[3] - x[2]) * (inputData[2] + 4 * inputData[3] + inputData[4]);
          halfBinComposite += lowerSimpson;
          break;
      }
    }
    return halfBinComposite;
  }

  /**
   * leftBoole.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double leftBoole(double[] inputData) {
    int n = inputData.length;
    double normalBoole = 0;
    for (int a = 1; a < n - 5; a += 4) {
      double area =
          (2.0 / 45.0)
              * (x[a + 1] - x[a])
              * (7 * inputData[a]
              + 32 * inputData[a + 1]
              + 12 * inputData[a + 2]
              + 32 * inputData[a + 3]
              + 7 * inputData[a + 4]);
      normalBoole += area;
    }

    return normalBoole;
  }

  /**
   * leftRectangularMethod.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double leftRectangularMethod(double[] inputData) {
    int n = inputData.length;
    double[] y = generateTestData_v1();
    double rectangularIntegral = 0;
    for (int a = 0; a < n - 2; a++) {
      double area = (x[a + 1] - x[a]) * y[a];
      rectangularIntegral += area;
    }
    return rectangularIntegral;
  }

  /**
   * leftSimpsons.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double leftSimpsons(double[] inputData) {
    int n = inputData.length;
    double normalSimpsons = 0;
    for (int a = 1; a < n - 4; a += 2) {
      double area =
          (1.0 / 3.0)
              * (x[a + 1] - x[a])
              * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
      normalSimpsons += area;
    }
    return normalSimpsons;
  }

  /**
   * leftTrapInput.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double leftTrapInput(double[] inputData) {
    int n = x.length;
    double sum = 0;
    double total = 0;
    double trapIntegral = 0;
    for (int a = 0; a < n - 2; a++) {
      if (a > 0) {
        double area = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
        sum = area + total;
        total = sum;
      }
      if (a == n - 3) {
        trapIntegral = sum;
      }
      if (a == 0) {
        total = (inputData[a + 1] + inputData[a]) / (double) 2 * (x[a + 1] - x[a]);
      }
    }

    return trapIntegral;
  }

  /**
   * main.
   *
   * @param args an array of {@link java.lang.String} objects.
   */
  public static void main(String[] args) {
    double testAnswer;
    double testTrap, testTrapRight, avgTrap, avgTrapError;
    double testSimp, testSimpRight, avgSimp, avgSimpError;
    double testBoole, testBooleRight, avgBoole, avgBooleError;
    double testRect, testRectRight, avgRect, avgRectError;
    DecimalFormat decimalFormat = new DecimalFormat("#.00");

    testAnswer = 2.98393938659796;
    System.out.print(" The test case answer is " + testAnswer + "\n\n");

    testRect = leftRectangularMethod(generateTestData_v1());
    testRectRight = rightRectangularMethod(generateTestData_v1());
    avgRect = (testRect + testRectRight) / 2.0;
    avgRectError = abs(testAnswer - avgRect) / testAnswer * 100;

    testTrap = leftTrapInput(generateTestData_v1());
    testTrapRight = rightTrapInput(generateTestData_v1());
    avgTrap = (testTrap + testTrapRight) / 2.0;
    avgTrapError = abs(avgTrap - testAnswer) / testAnswer * 100;

    testSimp =
        leftSimpsons(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 1, "left");
    testSimpRight =
        rightSimpsons(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 1, "right");
    avgSimp = (testSimp + testSimpRight) / 2.0;
    avgSimpError = abs(testAnswer - avgSimp) / testAnswer * 100;

    testBoole =
        leftBoole(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 2, "left");
    testBooleRight =
        rightBoole(generateTestData_v1()) + halfBinComposite(generateTestData_v1(), 2, "right");
    avgBoole = (testBoole + testBooleRight) / 2.0;
    avgBooleError = abs(testAnswer - avgBoole) / testAnswer * 100;

    // Average integrals
    System.out.print(" Average integrals \n\n");
    System.out.print(" Average rectangular " + avgRect + "\n");
    System.out.print(" Average trap " + avgTrap + "\n");
    System.out.print(" Average Simpsons " + avgSimp + "\n");
    System.out.print(" Average Boole " + avgBoole + "\n");

    // Average integral errors
    System.out.print("\n Average integral error \n\n");
    System.out.print(" Average rectangular error " + decimalFormat.format(avgRectError) + "%\n");
    System.out.print(" Average Trapezoidal error  " + decimalFormat.format(avgTrapError) + "%\n");
    System.out.print(" Average Simpsons error " + decimalFormat.format(avgSimpError) + "%\n");
    System.out.print(" Average Boole error " + decimalFormat.format(avgBooleError) + "%\n");
  }

  /**
   * rightBoole.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double rightBoole(double[] inputData) {
    int n = inputData.length;
    double normalBoole = 0;
    for (int a = 4; a < n - 5; a += 4) {
      // Simpsons and trapezoid + lower bin on left
      double area =
          (2.0 / 45.0)
              * (x[a + 1] - x[a])
              * (7 * inputData[a]
              + 32 * inputData[a + 1]
              + 12 * inputData[a + 2]
              + 32 * inputData[a + 3]
              + 7 * inputData[a + 4]);
      normalBoole += area;
    }

    return normalBoole;
  }

  /**
   * rightRectangularMethod.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double rightRectangularMethod(double[] inputData) {
    int n = inputData.length;
    double rectangularIntegral = 0;
    double[] y = generateTestData_v1();
    for (int a = 1; a < n - 1; a++) {
      double area = (x[a + 1] - x[a]) * y[a];
      rectangularIntegral += area;
    }
    return rectangularIntegral;
  }

  /**
   * rightSimpsons.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double rightSimpsons(double[] inputData) {
    int n = inputData.length;
    double normalSimpsons = 0;
    for (int a = 2; a < n - 3; a += 2) {
      // Extra trap on lower edge so right edge of rightmost bin aligns with the upper half bin
      double area =
          (1.0 / 3.0)
              * (x[a + 1] - x[a])
              * (inputData[a] + 4 * inputData[a + 1] + inputData[a + 2]);
      normalSimpsons += area;
    }
    return normalSimpsons;
  }

  /**
   * rightTrapInput.
   *
   * @param inputData an array of double values.
   * @return a double.
   */
  public static double rightTrapInput(double[] inputData) {
    int n = x.length;
    double sum = 0;
    double total = 0;
    double trapIntegral = 0;
    for (int a = 1; a < n - 1; a++) {
      if (a > 1) {
        double area = (inputData[a + 1] + inputData[a]) / 2.0 * (x[a + 1] - x[a]);
        sum = area + total;
        total = sum;
      }
      if (a == n - 2) {
        trapIntegral = sum;
      }
      if (a == 1) {
        total = (inputData[a + 1] + inputData[a]) / 2.0 * (x[a + 1] - x[a]);
      }
    }

    return trapIntegral;
  }
}
