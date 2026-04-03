#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#    Author: Malte Siemers, Freie Universität Berlin
#
#    If you use this software or anything it produces for work to be published,
#    please cite:
#
#    Malte Siemers, Michalis Lazaratos, Konstantina Karathanou,
#    Federico Guerra, Leonid Brown, and Ana-Nicoleta Bondar.
#    Bridge: A graph-based algorithm to analyze dynamic H-bond networks
#    in membrane proteins, Journal of Chemical Theory and Computation, 2019.
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning, module="MDAnalysis.*")

from . import helpfunctions as _hf
import MDAnalysis as _MDAnalysis
from MDAnalysis import transformations
import pickle as cPickle
import copy as _cp
from pathlib import Path

# import matplotlib
# matplotlib.use('TKAgg', warn=False)
import matplotlib.pyplot as _plt

import pdb


class BasicFunctionality(object):

    def __init__(
        self,
        selection=None,
        structure=None,
        trajectories=None,
        start=None,
        stop=None,
        step=1,
        ions=[],
        restore_filename=None,
        wrap_dcd=None,  # options are: None, 'fast', 'centered'
        write_wrapped_traj_to=None,
    ):

        if restore_filename is not None:
            self.load_from_file(restore_filename)
            return

        if selection is None:
            raise AssertionError("No selection string.")
        if structure is None:
            raise AssertionError("No structure file path.")
        self._selection = selection
        self._structure = structure
        self._trajectories = trajectories
        if trajectories is not None:
            self._universe = _MDAnalysis.Universe(structure, trajectories)
            if wrap_dcd == "fast":
                print("Wrapping trajectories with fast, simple method")
                self._universe.trajectory.add_transformations(
                    transformations.wrap(self._universe.atoms)
                )
            elif wrap_dcd == "centered":
                print("Wrapping trajectories with centered method")
                protein = self._universe.select_atoms("protein")
                solvent = self._universe.select_atoms("not protein")

                transforms = [
                    transformations.unwrap(self._universe.atoms),
                    transformations.center_in_box(protein, wrap=False),
                    transformations.wrap(solvent, compound="residues"),
                ]

                self._universe.trajectory.add_transformations(*transforms)

            if (wrap_dcd == "fast" or wrap_dcd == "centered") and write_wrapped_traj_to:
                lengths = [len(traj) for traj in self._universe.trajectory.readers]
                write_start = 0

                mode = "_centered" if wrap_dcd == "centered" else ""
                for i, length in enumerate(lengths):

                    input_path = Path(trajectories[i])

                    file_name = f"{input_path.stem}_PBC{mode}{input_path.suffix}"
                    output_path = Path(write_wrapped_traj_to, file_name)

                    write_stop = write_start + length
                    with _MDAnalysis.Writer(
                        str(output_path), n_atoms=self._universe.atoms.n_atoms
                    ) as W:
                        for ts in self._universe.trajectory[write_start:write_stop]:
                            W.write(self._universe.atoms)

                    write_start = write_stop

        else:
            self._universe = _MDAnalysis.Universe(structure)
        self._trajectory_slice = slice(
            start if isinstance(start, int) else None,
            stop if isinstance(stop, int) else None,
            step,
        )

        self._mda_selection = self._universe.select_atoms(selection)
        if not self._mda_selection:
            raise AssertionError("No atoms match the selection")

        self._ions_list = ions
        self._ions = _hf.EmptyGroup()
        self._ions_ids = []
        self._ions_ids_atomwise = []
        if ions:
            self._ions = self._universe.select_atoms(
                " or ".join(["name {}".format(ion) for ion in ions])
            )
            self._ions_ids = _hf.MDA_info_list(self._ions, detailed_info=False)
            self._ions_ids_atomwise = _hf.MDA_info_list(self._ions, detailed_info=True)

        self._water = self._universe.select_atoms(_hf.water_definition)
        self._water_ids = _hf.MDA_info_list(self._water, detailed_info=False)
        self._water_ids_atomwise = _hf.MDA_info_list(self._water, detailed_info=True)

        self.initial_results = {}
        self.filtered_results = {}
        self.nb_frames = len(
            [0 for i in self._universe.trajectory[self._trajectory_slice]]
        )
        self.applied_filters = {
            "resnames": None,
            "segnames": None,
            "shortest_paths": None,
            "single_path": None,
            "occupancy": None,
            "connected_component": None,
            "shells": None,
            "avg_least_bonds": None,
            "backbone": None,
        }

    def dump_to_file(self, fname):
        tmp = _cp.copy(self)
        tmp._universe = None
        tmp._water = None
        tmp._hydrogen = None
        tmp._mda_selection = None
        tmp._da_selection = None
        tmp._donors = None
        tmp._acceptors = None
        with open(fname, "wb") as af:
            af.write(cPickle.dumps(tmp.__dict__))

    def load_from_file(self, fname, reload_universe=False):
        with open(fname, "rb") as af:
            self.__dict__ = cPickle.loads(af.read())
        if reload_universe:
            self._reload_universe()

    def _reload_universe(self):
        if self._trajectories is not None:
            self._universe = _MDAnalysis.Universe(self._structure, self._trajectories)
        else:
            self._universe = _MDAnalysis.Universe(self._structure)
        self._mda_selection = self._universe.select_atoms(self._selection)
        self._water = self._universe.select_atoms(_hf.water_definition)

    def _set_results(self, result):
        self.initial_results = result
        self.filtered_results = result
        self._generate_graph_from_current_results()
        self._generate_filtered_graph_from_filtered_results()

    def _add_overwrite_results(self, result):
        for bond in result:
            self.initial_results[bond] = result[bond]
        self.filtered_results = self.initial_results
        self._generate_graph_from_current_results()
        self._generate_filtered_graph_from_filtered_results()

    def _save_or_draw(self, filename, data=None, return_figure=False):
        if filename is not None:
            end = filename.split(".")[-1]
            if end == "eps":
                _plt.text.usetex = True
            _plt.savefig(filename, format=end, dpi=300)
            _plt.close()
        else:
            _plt.show()

    def duplicate(self):
        return _cp.copy(self)
