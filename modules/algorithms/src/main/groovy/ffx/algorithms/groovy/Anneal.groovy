//******************************************************************************
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
//******************************************************************************
package ffx.algorithms.groovy

import ffx.algorithms.cli.AlgorithmsScript
import ffx.algorithms.cli.AnnealOptions
import ffx.algorithms.cli.DynamicsOptions
import ffx.algorithms.optimize.anneal.SimulatedAnnealing
import ffx.numerics.Potential
import ffx.potential.cli.AtomSelectionOptions
import ffx.potential.cli.WriteoutOptions
import org.apache.commons.io.FilenameUtils
import picocli.CommandLine.Command
import picocli.CommandLine.Mixin
import picocli.CommandLine.Parameters

/**
 * The Anneal script.
 * <br>
 * Usage:
 * <br>
 * ffxc Anneal [options] &lt;filename&gt;
 */
@Command(description = " Run simulated annealing on a system.", name = "Anneal")
class Anneal extends AlgorithmsScript {

  @Mixin
  AtomSelectionOptions atomSelectionOptions

  @Mixin
  DynamicsOptions dynamics

  @Mixin
  AnnealOptions anneal

  @Mixin
  WriteoutOptions writeOut

  /**
   * An XYZ or PDB input file.
   */
  @Parameters(arity = "1", paramLabel = "file",
      description = "XYZ or PDB input file.")
  private String filename

  private SimulatedAnnealing simulatedAnnealing = null

  private Potential potential

  /**
   * Anneal Constructor.
   */
  Anneal() {
    this(new Binding())
  }

  /**
   * Anneal Constructor.
   * @param binding The Groovy Binding to use.
   */
  Anneal(Binding binding) {
    super(binding)
  }

  /**
   * {@inheritDoc}
   */
  @Override
  Anneal run() {

    // Init the context and bind variables.
    if (!init()) {
      return this
    }

    // Init DynamicsOptions (e.g. the thermostat and barostat flags).
    dynamics.init()

    // Load the MolecularAssembly.
    activeAssembly = getActiveAssembly(filename)
    if (activeAssembly == null) {
      logger.info(helpString())
      return this
    }

    // Set the filename.
    filename = activeAssembly.getFile().getAbsolutePath()

    // Set active atoms.
    atomSelectionOptions.setActiveAtoms(activeAssembly)

    logger.info("\n Running simulated annealing on " + filename + "\n")

    // Define the path to the restart File.
    File dyn = new File(FilenameUtils.removeExtension(filename) + ".dyn")
    if (!dyn.exists()) {
      dyn = null
    }

    potential = activeAssembly.getPotentialEnergy()

    simulatedAnnealing = anneal.createAnnealer(dynamics, activeAssembly,
        activeAssembly.getPotentialEnergy(), algorithmListener, dyn)

    simulatedAnnealing.setPrintInterval(dynamics.report)
    simulatedAnnealing.setSaveFrequency(dynamics.write)
    simulatedAnnealing.setRestartFrequency(dynamics.checkpoint)
    simulatedAnnealing.setTrajectorySteps(dynamics.trajSteps)

    simulatedAnnealing.anneal()

    if (baseDir == null || !baseDir.exists() || !baseDir.isDirectory() || !baseDir.canWrite()) {
      baseDir = new File(FilenameUtils.getFullPath(filename))
    }

    String dirName = baseDir.toString() + File.separator
    String fileName = FilenameUtils.getName(filename)
    fileName = FilenameUtils.removeExtension(fileName)

    writeOut.saveFile(String.format("%s%s", dirName, fileName), algorithmFunctions, activeAssembly)

    return this
  }

  SimulatedAnnealing getAnnealing() {
    return simulatedAnnealing
  }

  /**
   * {@inheritDoc}
   */
  @Override
  List<Potential> getPotentials() {
    List<Potential> potentials
    if (potential == null) {
      potentials = Collections.emptyList()
    } else {
      potentials = Collections.singletonList(potential)
    }
    return potentials
  }

  /**
   * {@inheritDoc}
   */
  @Override
  boolean destroyPotentials() {
    return getPotentials().stream().allMatch({
      it.destroy()
    })
  }
}
