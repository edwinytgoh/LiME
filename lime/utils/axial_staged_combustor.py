import cantera as ct
import numpy as np

from .cantera_utils import correct_nox, mix, premix
from .flow_rates import solve_mass_flow_symbolic


class StagedCombustor(object):
    def __init__(self, args):
        """
        Initialize the StagedCombustor.
        
        Parameters:
        - args: argparse.Namespace object with the following attributes (specified in config.yaml):
            - phi_main: float, equivalence ratio of main burner
            - main_fuel: str, fuel name for main burner
            - oxidizer: str, oxidizer name for main and secondary burner
            - mech: str, path to mechanism file
            - P_atm: float, pressure of the combustor in atm
            - T_fuel_main: float, temperature of main fuel
            - T_oxidizer: float, temperature of main oxidizer
            - tau_secondary_ms: float, time at which to switch to secondary burner
            - phi_secondary: float, equivalence ratio of secondary burner
            - secondary_fuel: str, fuel name for secondary burner
            - T_fuel_secondary: float, temperature of secondary fuel
            - flame_width_m: float, width of the flame in meters
            - flame_refine_criteria: dict, refinement criteria for the flame object
            - flame_transport_model: str, transport model for the flame object
        """        
        self.parse_args(args)  # add each item in args as an attribute
        self.main_burner_gas = premix(
            args.phi_main,
            args.main_fuel,
            args.oxidizer,
            args.mech,
            args.P_atm * ct.one_atm,
            args.T_fuel_main,
            args.T_oxidizer,
        )
        self.main_stage_completed = False  # whether main stage has finished running
        self.sec_stage_completed = False  # whether secondary stage has finished

        self.initialize_flame(args)  # initialize flame object, but don't run yet

        self.column_names = (
            ["time", "x", "velocity", "T"]
            + ["Y_" + sn for sn in self.main_burner_gas.species_names]
            + ["X_" + sn for sn in self.main_burner_gas.species_names]
        )
        self.column_index = {col: i for i, col in enumerate(self.column_names)}

    def __getitem__(self, species):
        """
        Get the mole fraction of a species in the gas mixture at the current time.
        
        Parameters:
        - species: str, name of the species to get the mole fraction for
        
        Returns:
        - mole fraction of the species
        """        
        is_ppmvd = "_ppmvd" in species
        species = species.split("_ppmvd")[0]

        if species not in self.main_burner_gas.species_names:
            raise ValueError(f"Species {species} not in main burner gas")
        else:
            i = self.main_burner_gas.species_index(species)

        g = self.secondary_gas if self.sec_stage_completed else self.main_burner_gas

        if not is_ppmvd:
            return g.X[i]
        else:
            return correct_nox(self[species], self["H2O"], self["O2"])

    @property
    def T(self):
        """
        Get the temperature of the gas mixture at the current time.
        
        Returns:
        - temperature of the gas mixture
        """        
        return (
            self.secondary_gas.T if self.sec_stage_completed else self.main_burner_gas.T
        )

    @property
    def NO_ppmvd(self):
        """
        Get the concentration of NO in parts per million by volume dry at the current time.
        
        Returns:
        - concentration of NO in parts per million by volume dry
        """        
        return correct_nox(self["NO"], self["H2O"], self["O2"])

    def parse_args(self, args):
        """
        Add each item in args to self.__dict__, so that we can access them as attributes.
        
        Parameters:
        - args: argparse.Namespace object with the following attributes:
            - phi_main: float, equivalence ratio of main burner
            - main_fuel: str, fuel name for main burner
            - oxidizer: str, oxidizer name for main and secondary burner
            - mech: str, path to mechanism file
            - P_atm: float, pressure of the combustor in atm
            - T_fuel_main: float, temperature of main fuel
            - T_oxidizer: float, temperature of main oxidizer
            - tau_secondary_ms: float, time at which to switch to secondary burner
            - phi_secondary: float, equivalence ratio of secondary burner
            - secondary_fuel: str, fuel name for secondary burner
            - T_fuel_secondary: float, temperature of secondary fuel
        """
        for key, value in args.__dict__.items():
            self.__dict__[key] = value

        self.P = self.P_atm * ct.one_atm
        self.tau_secondary = self.tau_secondary_ms * 1e-3

        self.mfm, self.mom, self.mfs, self.mos = solve_mass_flow_symbolic(args)

    def initialize_flame(self, args):
        """Initialize the flame object.

        Parameters:
        - args: argparse.Namespace object with the following attributes:
            - flame_width_m: float, width of the flame in meters
            - flame_refine_criteria: dict, refinement criteria for the flame object
            - flame_transport_model: str, transport model for the flame object
        """
        self.flame = ct.FreeFlame(self.main_burner_gas, width=args.flame_width_m)
        self.flame_time = None
        self.flame.set_refine_criteria(**args.flame_refine_criteria)
        self.flame.transport_model = args.flame_transport_model

    def initialize_secondary_gas(self):
        """
        Initialize the gas mixture for the secondary burner using the specified temperature, pressure,
        fuel, and oxidizer.
        """        
        sec_fuel = ct.Solution(self.mech)
        sec_ox = ct.Solution(self.mech)
        sec_fuel.TPX = self.T_fuel_secondary, self.P, self.secondary_fuel
        sec_ox.TPX = self.T_oxidizer, self.P, self.oxidizer
        self.secondary_gas = mix(
            [sec_fuel, sec_ox], [self.mfs, self.mos], mech=self.mech, P=self.P
        )

    def run_main_burner(self):
        """Run the main burner and store the flame object and flame array. The flame array contains the
        time, axial position, velocity, and temperature of the flame, as well as the mole fractions
        of all species in the main burner gas mixture.
        """   
        self.flame_array = self.run_flame()
        vit_array = self.get_flame_state_at_time(self.tau_main_ms * 1e-3)

        # set main burner gas state to the last point in the vitiated gas array
        self.main_burner_gas.TPX = (
            vit_array[-1, self.column_index["T"]],  # Temperature
            self.flame.gas.P,  # Pressure
            vit_array[-1, slice(-self.flame.gas.n_species, None)],  # Mole fractions
        )

        self.main_burner_stage_completed = True
        return vit_array

    def run_flame(self):
        """Run the flame object and return the flame array. The flame array contains the time, axial
        position, velocity, and temperature of the flame, as well as the mole fractions of all species
        in the main burner gas mixture.
    
        Returns:
            np.ndarray: numpy array, with rows representing time points and columns representing the
            variables in self.column_names
        """
        # Run only if flame has not been run before,
        # i.e., if flame_time hasn't been set yet.
        if self.flame_time is None:
            f = self.flame
            f.solve(loglevel=0, auto=True, refine_grid=True)

            # To convert distance into time we need: 1) a t=0 point, and 2) a velocity
            # For t=0, we use the point of maximum CH2O. (or some other indicator species)
            ignition_species_index = f.gas.species_names.index(
                self.ignition_indicator_species
            )
            X_ignition = f.X[ignition_species_index]  # X_ign.shape = (n_grid, )
            start_idx = np.argmax(X_ignition)
            average_velocity = (
                np.array(f.velocity[start_idx:] + f.velocity[start_idx - 1 : -1]) * 0.5
            )
            dx = np.array(f.grid[start_idx:] - f.grid[start_idx - 1 : -1])
            dt = dx / average_velocity

            # add pre-time to "fill out" time array before start of the flame
            # (because distance = 0 corresponds to negative time, i.e., time "before" the flame)
            pre_time = [-999] * (start_idx - 1)
            pre_time.extend([0])
            # numerically integrate (read: sum) dt to get array of times for each x location
            self.flame_time = np.hstack((np.array(pre_time), np.cumsum(dt)))

        # Return flame state at each time step
        flame_vars = [
            [self.flame_time, self.flame.grid, self.flame.velocity, self.flame.T],
            self.flame.Y,
            self.flame.X,
        ]
        self.flame_array = np.concatenate(flame_vars, axis=0).T
        return self.flame_array

    def get_flame_state_at_time(self, tau_required_s):
        """
        Returns the state of the flame at a given time - `tau_required`.
        For a given flame_array, there are 2 cases that we need to consider:
        1. `tau_required` < max flame time:
            Find the closest time in flame_array and decide whether to interpolate
        2. `tau_required` > max flame time:
            Start the vitiator at the end of the flame simulation and run it for the
            remaining time.

        For case 1, interpolate (i.e., run vitiator) if the time difference between the closest time
        and `tau_required` is less than a certain threshold.

        """
        time_col_idx = self.column_index["time"]
        temp_col_idx = self.column_index["T"]
        velocity_col_idx = self.column_index["velocity"]
        X_col_idx = slice(-self.flame.gas.n_species, None)

        flame_time = self.flame_array[:, time_col_idx]
        flame_temp = self.flame_array[:, temp_col_idx]
        flame_X = self.flame_array[:, X_col_idx]

        flame_end_idx, tau_vit = self.get_tau_vit(tau_required_s, flame_time)
        vit_array = self.flame_array[0 : flame_end_idx + 1, :]

        if tau_vit > 0:
            # Create gas and run batch reactor to "extend" the flame
            v_gas = ct.Solution(self.mech)
            v_gas.TPX = (
                flame_temp[flame_end_idx],
                self.flame.gas.P,
                flame_X[flame_end_idx, :],
            )
            vitiator = ct.ConstPressureReactor(v_gas)
            vit_network = ct.ReactorNet([vitiator])
            Δt = 1e-7  # time step for interpolation
            for t in np.arange(0, tau_vit, Δt):
                vit_network.advance(t)

            # add data at interpolated point to NumPy array
            t_ = t + flame_time[flame_end_idx]  # time after vitiator
            x_ = self.flame.grid[flame_end_idx]  # dummy distance after flame
            u_ = vit_array[-1, velocity_col_idx]  # dummy velocity after flame
            new_row = np.concatenate([[t_, x_, u_, v_gas.T], v_gas.Y, v_gas.X])
            vit_array = np.concatenate((vit_array, new_row[np.newaxis, :]), axis=0)
        self.vit_array = vit_array

        return vit_array
        # DataFrame with final data (useful for sending to csv later)
        # column_names = ["t", "T"] + ["X_" + sn for sn in species_names]
        # flame_data_frame = pd.DataFrame(data=vit_array.transpose(), columns=column_names)
        # return vitiator_gas, flame_data_frame

    def get_tau_vit(self, tau_required, flame_time, interp_threshold=1e-5):
        """
        Determine the index of the flame time array from which to attach a vitiator.

        Args:
            tau_required (float): required flame end time in seconds
            flame_time (np.array): array of flame times in seconds
            interp_threshold (float): threshold to decide whether to run vitiator or not. Default is 1e-5.

        Returns:
            vitiator_start_index (int): index of flame time array from which to attach vitiator
            tau_vit (float): vitiator duration in seconds. If 0 then no vitiator is required.
        """
        # if tau_required >= current residence time, figure out how long to run the vitiator for (tau_vit)
        if tau_required >= max(flame_time):
            vitiator_start_index = -1  # start from the end of the flame simulation
            tau_vit = tau_required - max(flame_time)
        else:
            # find index at which tau_required is closest to an entry in t:
            vitiator_start_index = np.abs(flame_time - tau_required).argmin()
            # Run vitiator if closest time differs from required tau by more than 10 microseconds
            if abs(flame_time[vitiator_start_index] - tau_required) > interp_threshold:
                # need to interpolate from a point lower than tau_required
                vitiator_start_index -= 1
                tau_vit = tau_required - flame_time[vitiator_start_index]
            else:  # don't need to vitiate if we're already very close to the required time
                tau_vit = 0
        return vitiator_start_index, tau_vit

    def save_main_burner_to_file(self, filename):
        # main_burner_df = pd.DataFrame(
        #     data=flame_array, columns=["Time"] + self.column_names,
        # )
        # main_burner_df["P"] = self.flame.P
        # # Write dataframe to file in parquet format
        # table = pa.Table.from_pandas(main_burner_df)
        # main_burner_out_file = os.path.join(
        #     self.output_path, "{0}_{1:.4f}.pickle".format("phi_main", self.phi_main)
        # )
        # pq.write_table(table, main_burner_out_file)
        pass

    def run_secondary_stage(self):
        """
        Run the secondary burner and store the secondary gas mixture and secondary flame array. The
        secondary flame array contains the time, axial position, velocity, and temperature of the
        flame, as well as the mole fractions of all species in the secondary gas mixture.
        
        Parameters:
        - None
        
        Returns:
        - None
        """        
        # Mix main burner gas with secondary reactants
        sec_fuel = ct.Solution(self.mech)
        sec_ox = ct.Solution(self.mech)
        sec_fuel.TPX = self.T_fuel_secondary, self.P, self.secondary_fuel
        sec_ox.TPX = self.T_oxidizer, self.P, self.oxidizer
        self.secondary_gas = mix(
            [sec_fuel, sec_ox, self.main_burner_gas],
            [self.mfs, self.mos, self.mfm + self.mom],
            mech=self.mech,
            P=self.P,
        )
        self.secondary_gas.name = "Secondary gas"
        print(self.secondary_gas())
        self.secondary_stage = ct.ConstPressureReactor(self.secondary_gas)
        network = ct.ReactorNet([self.secondary_stage])
        time_array = np.arange(0, self.tau_secondary, self.secondary_dt)
        for t in time_array:
            network.advance(t)
        print(self.secondary_gas())
        self.sec_stage_completed = True
