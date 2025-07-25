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

import ffx.algorithms.AlgorithmListener;
import ffx.algorithms.dynamics.MolecularDynamics;
import ffx.crystal.CrystalPotential;
import ffx.potential.MolecularAssembly;
import ffx.potential.bonded.LambdaInterface;
import ffx.potential.cli.WriteoutOptions;
import picocli.CommandLine.ArgGroup;
import picocli.CommandLine.Option;

import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Logger;

import static java.lang.String.format;

/**
 * Represents command line options for scripts that calculate thermodynamics.
 *
 * @author Michael J. Schnieders
 * @author Jacob M. Litman
 * @since 1.0
 */
public class ThermodynamicsOptions {

  private static final Logger logger = Logger.getLogger(ThermodynamicsOptions.class.getName());

  /**
   * The ArgGroup keeps the ThermodynamicsOptions together when printing help.
   */
  @ArgGroup(heading = "%n Thermodynamics Options%n", validate = false)
  private final ThermodynamicsOptionGroup group = new ThermodynamicsOptionGroup();

  /**
   * Return the selected Thermodynamics algorithm as an enumerated type.
   *
   * @return Corresponding thermodynamics algorithm
   */
  public ThermodynamicsAlgorithm getAlgorithm() {
    return ThermodynamicsAlgorithm.parse(group.thermoAlgoString);
  }

  /**
   * getEquilSteps.
   *
   * @return The number of equilibration steps.
   */
  public long getEquilSteps() {
    return group.equilibrationSteps;
  }

  /**
   * Getter for the field <code>resetNumSteps</code>.
   *
   * @return a boolean.
   */
  public boolean getResetNumSteps() {
    return group.resetNumSteps;
  }

  /**
   * Run an alchemical free energy window.
   *
   * @param molecularAssemblies All involved MolecularAssemblies.
   * @param crystalPotential    The Potential to be sampled.
   * @param dynamicsOptions     DynamicsOptions.
   * @param writeoutOptions     WriteoutOptions
   * @param dyn                 MD restart file
   * @param algorithmListener   AlgorithmListener
   * @return The MolecularDynamics object constructed.
   */
  public MolecularDynamics runFixedAlchemy(MolecularAssembly[] molecularAssemblies,
                                           CrystalPotential crystalPotential, DynamicsOptions dynamicsOptions,
                                           WriteoutOptions writeoutOptions, File dyn, AlgorithmListener algorithmListener) {
    dynamicsOptions.init();

    MolecularDynamics molDyn = dynamicsOptions.getDynamics(writeoutOptions, crystalPotential,
        molecularAssemblies[0], algorithmListener);
    for (int i = 1; i < molecularAssemblies.length; i++) {
      molDyn.addAssembly(molecularAssemblies[i]);
    }

    boolean initVelocities = true;
    long nSteps = dynamicsOptions.getSteps();
    molDyn.setRestartFrequency(dynamicsOptions.getCheckpoint());
    // Start sampling.
    if (group.equilibrationSteps > 0) {
      logger.info("\n Beginning Equilibration");
      runDynamics(molDyn, group.equilibrationSteps, dynamicsOptions, writeoutOptions, true, dyn);
      logger.info(" Beginning Fixed-Lambda Alchemical Sampling");
      initVelocities = false;
    } else {
      logger.info("\n Beginning Fixed-Lambda Alchemical Sampling Without Equilibration");
      if (!group.resetNumSteps) {
        // Workaround for being unable to pick up pre-existing steps.
        initVelocities = true;
      }
    }

    if (nSteps > 0) {
      runDynamics(molDyn, nSteps, dynamicsOptions, writeoutOptions, initVelocities, dyn);
    } else {
      logger.info(" No steps remaining for this process!");
    }
    return molDyn;
  }

  /**
   * Run a non-equilibrium alchemical free energy simulation.
   *
   * @param molecularAssemblies All involved MolecularAssemblies.
   * @param crystalPotential    The Potential to be sampled.
   * @param dynamicsOptions     DynamicsOptions.
   * @param writeoutOptions     WriteoutOptions
   * @param dyn                 MD restart file
   * @param algorithmListener   AlgorithmListener
   * @return The MolecularDynamics object constructed.
   */
  public MolecularDynamics runNEQ(MolecularAssembly[] molecularAssemblies,
                                  CrystalPotential crystalPotential, DynamicsOptions dynamicsOptions,
                                  WriteoutOptions writeoutOptions, File dyn, AlgorithmListener algorithmListener) {
    dynamicsOptions.init();

    MolecularDynamics molDyn = dynamicsOptions.getDynamics(writeoutOptions, crystalPotential,
        molecularAssemblies[0], algorithmListener);
    for (int i = 1; i < molecularAssemblies.length; i++) {
      molDyn.addAssembly(molecularAssemblies[i]);
    }

    boolean initVelocities = true;
    long nSteps = dynamicsOptions.getSteps();
    molDyn.setRestartFrequency(dynamicsOptions.getCheckpoint());
    // Start sampling.
    if (group.equilibrationSteps > 0) {
      double initialLambda = group.reverseNEQ ? 1.0 : 0.0;
      logger.info(format("\n Beginning Equilibration (at L=%5.3f)", initialLambda));
      LambdaInterface lambdaInterface = (LambdaInterface) crystalPotential;
      lambdaInterface.setLambda(initialLambda);
      runDynamics(molDyn, group.equilibrationSteps, dynamicsOptions, writeoutOptions, true, dyn);
      if (nSteps > 0) {
        logger.info(" Beginning Non-Equilibrium Sampling");
      }
      initVelocities = false;
    } else if (nSteps > 0) {
      logger.info("\n Beginning Non-Equilibrium Sampling Without Equilibration");
      if (dyn.exists()) {
        initVelocities = false;
      }
    }

    if (nSteps > 0) {
      molDyn.setNonEquilibriumLambda(true, group.nonEquilibriumSteps, group.reverseNEQ);
      runDynamics(molDyn, nSteps, dynamicsOptions, writeoutOptions, initVelocities, dyn);
    }

    return molDyn;
  }

  /**
   * The number of equilibration steps prior to production OST counts begin.
   *
   * @return Returns the number of equilibration steps.
   */
  public long getEquilibrationSteps() {
    return group.equilibrationSteps;
  }

  private void runDynamics(MolecularDynamics molecularDynamics, long nSteps,
                           DynamicsOptions dynamicsOptions, WriteoutOptions writeoutOptions, boolean initVelocities,
                           File dyn) {
    molecularDynamics.dynamic(nSteps, dynamicsOptions.getDt(), dynamicsOptions.getReport(),
        dynamicsOptions.getWrite(), dynamicsOptions.getTemperature(), initVelocities,
        writeoutOptions.getFileType(), dynamicsOptions.getCheckpoint(), dyn);
  }

  public void setEquilibrationSteps(long equilibrationSteps) {
    group.equilibrationSteps = equilibrationSteps;
  }

  /**
   * Ignores steps detected in .lam lambda-restart files.
   *
   * @return Returns true if the number of steps is being reset.
   */
  public boolean isResetNumSteps() {
    return group.resetNumSteps;
  }

  public void setResetNumSteps(boolean resetNumSteps) {
    group.resetNumSteps = resetNumSteps;
  }

  /**
   * The algorithm to be used (e.g. OST, window-based methods, etc.).
   *
   * @return Returns a String for requested algorithm.
   */
  public String getThermoAlgoString() {
    return group.thermoAlgoString;
  }

  public void setThermoAlgoString(String thermoAlgoString) {
    group.thermoAlgoString = thermoAlgoString;
  }

  /**
   * Collection of Thermodynamics Options.
   */
  private static class ThermodynamicsOptionGroup {

    /**
     * -Q or --equilibrate sets the number of equilibration steps prior to production OST counts
     * begin.
     */
    @Option(names = {"-Q", "--equilibrate"}, paramLabel = "1000", defaultValue = "1000",
        description = "Number of equilibration steps before evaluation of thermodynamics.")
    private long equilibrationSteps = 1000;

    /**
     * --nEQ or --nonEquilibriumSteps sets the number of non-equilibrium lambda steps.
     */
    @Option(names = {"--nEQ", "--nonEquilibriumSteps"}, paramLabel = "100", defaultValue = "100",
        description = "Sets the number of non-equilibrium lambda steps.")
    private int nonEquilibriumSteps = 100;

    /**
     * --rNEQ or --reverseNonEquilibrium Run non-equilibrium dynamics from L=1 to L=0.
     */
    @Option(names = {"--rNEQ", "--reverseNonEquilibrium"}, defaultValue = "false",
        description = "Run non-equilibrium dynamics from L=1 to L=0.")
    private boolean reverseNEQ = false;

    /**
     * -rn or --resetNumSteps, ignores steps detected in .lam lambda-restart files and thus resets
     * the histogram; use -rn false to continue from the end of any prior simulation.
     */
    @Option(names = {"--rn", "--resetNumSteps"}, defaultValue = "false",
        description = "Ignore prior steps logged in .lam or similar files.")
    private boolean resetNumSteps = false;

    /**
     * --tA or --thermodynamicsAlgorithm specifies the algorithm to be used; currently serves as a
     * switch between OST, Fixed and NEQ methods.
     */
    @Option(names = {"--tA", "--thermodynamicsAlgorithm"}, paramLabel = "OST", defaultValue = "OST",
        description = "Choice of thermodynamics algorithm [OST, FIXED, or NEQ].")
    private String thermoAlgoString = "OST";
  }

  /**
   * Represents categories of thermodynamics algorithms that must be handled differently.
   */
  public enum ThermodynamicsAlgorithm {
    OST("OST", "MC-OST", "MD-OST", "DEFAULT"),
    FIXED("FIXED", "BAR", "MBAR", "FEP", "WINDOWED"),
    NEQ("NEQ", "NON-EQUILIBRIUM");

    private final Set<String> aliases;

    ThermodynamicsAlgorithm(String... aliases) {
      // If, for some reason, there are 100+ aliases, might change to a HashSet.
      Set<String> names = new TreeSet<>(Arrays.asList(aliases));
      names.add(this.name());
      this.aliases = Collections.unmodifiableSet(names);
    }

    /**
     * Parse a String to a corresponding thermodynamics algorithm, recognizing aliases.
     *
     * @param name Name to parse
     * @return A ThermodynamicsAlgorithm.
     * @throws IllegalArgumentException If name did not correspond to any alias of any
     *                                  ThermodynamicsAlgorithm.
     */
    public static ThermodynamicsAlgorithm parse(String name) throws IllegalArgumentException {
      String ucName = name.toUpperCase();
      for (ThermodynamicsAlgorithm thermodynamicsAlgorithm : values()) {
        if (thermodynamicsAlgorithm.aliases.contains(ucName)) {
          return thermodynamicsAlgorithm;
        }
      }
      throw new IllegalArgumentException(
          format(" Could not parse %s as a ThermodynamicsAlgorithm", name));
    }
  }
}
