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
package ffx.algorithms.cli;

import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

/**
 * Represents command line options for scripts that utilize a mobile lambda particle, such as
 * Thermodynamics.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class LambdaParticleOptions {

  /**
   * The default theta friction is 91 per psec.
   */
  private static final double DEFAULT_THETA_FRICTION = 91.0;

  /**
   * The default theta mass is 10 AMU.
   */
  private static final double DEFAULT_THETA_MASS = 10.0;

  /**
   * The ArgGroup keeps the LambdaParticleOptions together when printing help.
   */
  @ArgGroup(heading = "%n Lambda Particle Options for MD-OST%n", validate = false)
  private final LambdaParticleOptionGroup group = new LambdaParticleOptionGroup();

  /**
   * Friction on the lambda particle.
   *
   * @return a double.
   */
  public double getLambdaFriction() {
    return group.lambdaFriction;
  }

  public void setLambdaFriction(double lambdaFriction) {
    group.lambdaFriction = lambdaFriction;
  }

  /**
   * The mass of the lambda particle.
   *
   * @return a double.
   */
  public double getLambdaMass() {
    return group.lambdaMass;
  }

  public void setLambdaMass(double lambdaMass) {
    group.lambdaMass = lambdaMass;
  }

  /**
   * Collection of Lambda Particle Options.
   */
  private static class LambdaParticleOptionGroup {

    /**
     * -m or --lambdaMass to set the mass of the lambda particle.
     */
    @Option(names = {"--lm", "--lambdaMass"}, paramLabel = "" + DEFAULT_THETA_MASS, defaultValue = "" + DEFAULT_THETA_MASS,
        description = "Mass of the lambda particle (AMU).")
    private double lambdaMass = DEFAULT_THETA_MASS;

    /**
     * -x or --lambdaFriction to set the friction on the lambda particle.
     */
    @Option(names = {"--lf", "--lambdaFriction"}, paramLabel = "" + DEFAULT_THETA_FRICTION, defaultValue = "" + DEFAULT_THETA_FRICTION,
        description = "Friction on the lambda particle (1/psec).")
    private double lambdaFriction = DEFAULT_THETA_FRICTION;
  }
}
