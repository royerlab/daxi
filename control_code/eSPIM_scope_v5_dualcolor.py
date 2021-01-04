"""making ao channel to send a sequence of signals"""
import nidaqmx
import nidaqmx.system
import numpy as np
import time


def set_dio_state(ch, value):
    """set a DIO channel,
    false to low, true to high"""
    with nidaqmx.Task() as task:
        task.do_channels.add_do_chan(ch)
        task.write([value], auto_start=True)


def set_ao_value(ch, value):
    """set the ao channel to value"""
    with nidaqmx.Task() as task:
        task.ao_channels.add_ao_voltage_chan(ch)
        task.write([value], auto_start=True)


class NIDaq:
    # Channel information
    ch_ao0 = "cDAQ1AO/ao0"  # scanning channel
    ch_ao1 = "cDAQ1AO/ao1"  # view switching
    ch_ao2 = "cDAQ1AO/ao2"  # view switching
    ch_ao3 = "cDAQ1AO/ao3"  # stripe reduction

    ch_ctr0 = "cDAQ1/_ctr0"  # for retrigger
    ch_ctr0_internal_output = "/cDAQ1/Ctr0InternalOutput"
    ch_ctr1 = "cDAQ1/_ctr1"  # for counting number of frames
    ch_ctr1_internal_output = "/cDAQ1/Ctr1InternalOutput"
    ch_ctr2 = "cDAQ1/_ctr2"  # idle
    ch_ctr3 = "cDAQ1/_ctr3"  # idle

    PFI0 = "/cDAQ1/PFI0"  # camera exposure input
    PFI1 = "/cDAQ1/PFI1"    # idle

    ch_dio0 = "cDAQ1DIO/port0/line0"  # 488 digital channel
    ch_dio1 = "cDAQ1DIO/port0/line1"  # 561 digital channel
    ch_dio2 = "cDAQ1DIO/port0/line2"  # bright field
    ch_dio3 = "cDAQ1DIO/port0/line3"  # idle
    ch_dio4 = "cDAQ1DIO/port0/line4"  # idle
    ch_dio5 = "cDAQ1DIO/port0/line5"  # idle
    ch_dio6 = "cDAQ1DIO/port0/line6"  # idle
    ch_dio7 = "cDAQ1DIO/port0/line7"  # idle

    # constants
    MAX_VOL = 10  # unit: v, maximum voltage of the ao channels
    MIN_VOL = -10  # unit: v, minimal voltage of the ao channels
    CONVERT_RATIO = 159  # unit: um / v, to convert from voltage to the scan distance of the galvo
    READOUT_TIME_FULL_CHIP = 0.01  # unit: second, the readout time of the full chip camera
    MAX_VERTICAL_PIXELS = 2048     # unit: pixels, maximal number of pixels along the vertical direction

    def __init__(
            self,
            exposure: float,
            nb_timepoints: int,
            scan_step: float,
            stage_scan_range: float,
            *,
            vertical_pixels: int = 2048,    # unit: pixels, number of pixels along the vertical direction
            num_samples: int = 20,  # for timing
            offset_view1: float  = 1550,   # unit: um, offset needed for view1
            offset_view2: float  = 1650,   # unit: um, offset needed for view2
            view1_galvo1: float  = 4.2,    # unit: v, to apply on galvo 1 for view 1
            view1_galvo2: float  = -4.08,  # unit: v, to apply on galvo 2 for view 1
            view2_galvo1: float  = -4.37,  # unit: v, to apply on galvo 1 for view 2
            view2_galvo2: float  = 3.66,   # unit: v, to apply on galvo 2 for view 2
            stripe_reduction_range: float = 0.1,    # unit: v, to apply on glavo gamma to reduce stripe
            stripe_reduction_offset: float = 0.58  # unit: v, to apply on glavo gamma to reduce stripe
        ):
        """
        Constructor of NIDaq.

        Parameters
        ----------
        exposure : float
            Expose time. Its unit is seconds.
        nb_timepoints : int
            Number of timepoints. Number of stacks to acquire.
        scan_step : float
            Scanning step size. unit: um.
        stage_scan_range : float
            Total range of stage scanning. unit: um.
        readout_time : float
             Readout time of the camera for full chip, unit: second.

        # TODO: add missing docstrings for rest of the optional arguments.
        """
        self.exposure = exposure
        self.nb_timepoints = nb_timepoints
        self.scan_step = scan_step
        self.stage_scan_range = stage_scan_range
        self.readout_time = self.READOUT_TIME_FULL_CHIP * vertical_pixels / self.MAX_VERTICAL_PIXELS
        self.num_samples = num_samples
        self.offset_view1 = offset_view1
        self.offset_view2 = offset_view2
        self.view1_galvo1 = view1_galvo1
        self.view1_galvo2 = view1_galvo2
        self.view2_galvo1 = view2_galvo1
        self.view2_galvo2 = view2_galvo2
        self.stripe_reduction_range = stripe_reduction_range
        self.stripe_reduction_offset = stripe_reduction_offset

        print("number of slices:", self.nb_slices)
        print("current sampling_rate is:", self.sampling_rate)

    @property
    def nb_slices(self):
        return int(self.stage_scan_range / self.scan_step) # for scanning, plus 2 due to Flash4 camera

    @property
    def sampling_rate(self):
        return self.num_samples / self.exposure

    def _get_ao_data(self, view: str):
        """generate the ndarray for an ao channel"""

        # for stripe reduction
        stripe_min = -self.stripe_reduction_range + self.stripe_reduction_offset
        stripe_max = self.stripe_reduction_range + self.stripe_reduction_offset
        nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
        nb_off_sample = round(self.readout_time * self.sampling_rate)
        data_ao3 = list(np.linspace(stripe_min, stripe_max, nb_on_sample))
        data_ao3.extend([stripe_min] * nb_off_sample)

        # for view switching and light sheet stabilization
        if view == "view1":
            offset = self._offset_dis_to_vol(self.offset_view1)   # convert the offset from um to v
            min_range = - self.scan_step / 2 / self.CONVERT_RATIO
            max_range = self.scan_step / 2 / self.CONVERT_RATIO
            data_ao0 = list(np.linspace(max_range + offset, min_range + offset, self.num_samples))
            data_ao1 = [self.view1_galvo1] * self.num_samples
            data_ao2 = [self.view1_galvo2] * self.num_samples
            return [data_ao0, data_ao1, data_ao2, data_ao3]
        elif view == "view2":
            offset = self._offset_dis_to_vol(self.offset_view2)  # convert the offset from um to v
            min_range = - self.scan_step / 2 / self.CONVERT_RATIO
            max_range = self.scan_step / 2 / self.CONVERT_RATIO
            data_ao0 = list(np.linspace(max_range + offset, min_range + offset, self.num_samples))
            data_ao1 = [self.view2_galvo1] * self.num_samples
            data_ao2 = [self.view2_galvo2] * self.num_samples
            return [data_ao0, data_ao1, data_ao2, data_ao3]

    def _get_do_data(self, channels):
        """
        Method to get digital output data.
        """
        if len(channels) == 1:
            nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
            data = [True] * nb_on_sample + [False] * (self.num_samples - nb_on_sample)
            return [data]
        elif len(channels) == 2:
            nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
            nb_off_sample = round(self.readout_time * self.sampling_rate)
            data_on = [True] * nb_on_sample + [False] * (self.num_samples - nb_on_sample)
            data_off = [False] * self.num_samples
            # return [[data_on, data_off], [data_off, data_on]]
            return [data_on + data_off, data_off + data_on]
        else:
            raise ValueError('Only supported up to 2 channels for now')

    def _crate_ao_task_for_acquisition(self):
        """create ao task for acquisition, include all the needed ao channels"""
        task_ao = nidaqmx.Task("ao0")
        task_ao.ao_channels.add_ao_voltage_chan(self.ch_ao0)
        task_ao.ao_channels.add_ao_voltage_chan(self.ch_ao1)
        task_ao.ao_channels.add_ao_voltage_chan(self.ch_ao2)
        task_ao.ao_channels.add_ao_voltage_chan(self.ch_ao3)
        return task_ao

    def acquire_stacks(self, channels, view):
        """acquire stackes, depending on the given channel and view.
        view=0, first view only;
        view=1, second view only;
        view=2, both views;
        channel=[488], [561], or [488, 561]
        """
        task_ao = self._crate_ao_task_for_acquisition()
        task_do = nidaqmx.Task("do0")  # for laser control

        # slect channel
        if channels == [488]:
            task_do.do_channels.add_do_chan(self.ch_dio0)
        elif channels == [561]:
            task_do.do_channels.add_do_chan(self.ch_dio1)
        elif channels == [488, 561]:
            task_do.do_channels.add_do_chan(self.ch_dio0)
            task_do.do_channels.add_do_chan(self.ch_dio1)
        else:
            raise ValueError('Channel not supported')

        # slect view
        if view == 0:
            self._acquire_stacks(task_ao, task_do, channels, ["view1"])
        elif view == 1:
            self._acquire_stacks(task_ao, task_do, channels, ["view2"])
        elif view == 2:
            self._acquire_stacks(task_ao, task_do, channels, ["view1", "view2"])
        else:
            raise ValueError('View not supported')

        task_ao.close()
        task_do.close()

    def _acquire_stacks(self, task_ao, task_do, channels, views):
        """set up the workflow to acquire multiple stacks"""
        data_ao = [self._get_ao_data(v) for v in views]  # different for each view due to different offsets
        data_do = self._get_do_data(channels)    # get the digital output data depending on the channels

        # set up the counter for loop through a zstack
        task_ctr_loop = nidaqmx.Task("counter0")
        ctr_loop = task_ctr_loop.ci_channels.add_ci_count_edges_chan(self.ch_ctr1, edge=nidaqmx.constants.Edge.RISING)
        ctr_loop.ci_count_edges_term = self.PFI0

        # set up the counter to retrigger the ao and do channels
        task_ctr_retrig = nidaqmx.Task("counter1")
        task_ctr_retrig.co_channels.add_co_pulse_chan_freq(self.ch_ctr0,
                                                           idle_state=nidaqmx.constants.Level.LOW,
                                                           freq=self.sampling_rate)
        task_ctr_retrig.timing.cfg_implicit_timing(sample_mode=nidaqmx.constants.AcquisitionType.FINITE,
                                                   samps_per_chan=self.num_samples)
        task_ctr_retrig.triggers.start_trigger.cfg_dig_edge_start_trig(trigger_source=self.PFI0,
                                                                       trigger_edge=nidaqmx.constants.Slope.RISING)
        task_ctr_retrig.triggers.start_trigger.retriggerable = True

        # set up ao channel
        task_ao.timing.cfg_samp_clk_timing(rate=self.sampling_rate,
                                           source=self.ch_ctr0_internal_output,
                                           sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)

        task_ctr_retrig.start()
        for _ in range(self.nb_timepoints):
            for v in range(len(data_ao)):  # change view
                # print("ao data size: ", len(data_ao))
                # print("ao data size2: ", len(data_ao[v]))
                # print("ao data size3: ", len(data_ao[v][0]))
                task_ao.write(data_ao[v])
                task_ao.start()
                for i_ch in range(1): #range(len(channels)):  # change channel # run once
                    # set up the do channel
                    task_do.timing.cfg_samp_clk_timing(rate=self.sampling_rate,
                                                       source=self.ch_ctr0_internal_output,
                                                       sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)
                    task_do.write(data_do)
                    task_do.start()

                    task_ctr_loop.start()
                    counts = 0
                    while counts < self.nb_slices + 1:  # Flash4.0 outputs 1 more pulse than asked
                        counts = task_ctr_loop.read()
                        time.sleep(0.005)  # wait time during loop
                    task_ctr_loop.stop()
                    print("counts: ", counts)
                    # print("counts: ", counts)
                    time.sleep(self.exposure + 0.1)  # add time to allow ao and do output for the last frame
                    task_do.stop()
                    print("one stack done!")
                task_ao.stop()

        # task_ctr_retrig.stop()
        task_ctr_retrig.close()
        task_ctr_loop.close()

    def _select_view(self, offset, galvo1, galvo2):
        """
        select the views, by using the correct offset of the scanning gavlo and the voltages of the switching galvos
        float -> void
        """
        offset_in_vol = self._offset_dis_to_vol(offset)
        data = [offset_in_vol, galvo1, galvo2]
        with nidaqmx.Task() as task:
            task.ao_channels.add_ao_voltage_chan(self.ch_ao0)
            task.ao_channels.add_ao_voltage_chan(self.ch_ao1)
            task.ao_channels.add_ao_voltage_chan(self.ch_ao2)
            task.write(data, auto_start=True)

    def select_view(self, view):
        """select one view in live mode"""
        if view == 1:
            self._select_view(self.offset_view1, self.view1_galvo1, self.view1_galvo2)
        elif view == 2:
            self._select_view(self.offset_view2, self.view2_galvo1, self.view2_galvo2)
        else:
            raise ValueError('View not supported')

    def _offset_dis_to_vol(self, offset):
        """
        convert the offset of each view from um to voltage
        float -> float
        um    -> v
        """
        return offset / self.CONVERT_RATIO - self.MAX_VOL

    def _select_channel(self, ch):
        nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
        # nb_off_sample = round(self.readout_time * self.sampling_rate)
        data = [True] * nb_on_sample + [False] * (self.num_samples - nb_on_sample)

        task = nidaqmx.Task()
        task.do_channels.add_do_chan(ch)
        self._retriggable_task(task, data)
        set_dio_state(ch, False)

    def _select_channel_remove_stripes(self, ch):
        nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
        # nb_off_sample = round(self.readout_time * self.sampling_rate)
        data_do = [True] * nb_on_sample + [False] * (self.num_samples - nb_on_sample)

        task_do = nidaqmx.Task()
        task_do.do_channels.add_do_chan(ch)

        task_ao = nidaqmx.Task("ao0")
        task_ao.ao_channels.add_ao_voltage_chan(self.ch_ao3)
        # for stripe reduction
        stripe_min = -self.stripe_reduction_range + self.stripe_reduction_offset
        stripe_max = self.stripe_reduction_range + self.stripe_reduction_offset
        nb_on_sample = round((self.exposure - self.readout_time) * self.sampling_rate)
        nb_off_sample = round(self.readout_time * self.sampling_rate)
        data_ao3 = list(np.linspace(stripe_min, stripe_max, nb_on_sample))
        data_ao3.extend([stripe_min] * nb_off_sample)

        self._retriggable_task_remove_stripes(task_do, data_do, task_ao, data_ao3)
        set_dio_state(ch, False)

    def select_channel(self, ch):
        """select one channel in live mode"""
        if ch == 488:
            self._select_channel(self.ch_dio0)
        elif ch == 561:
            self._select_channel(self.ch_dio1)
        else:
            raise ValueError('Channel not supported!')

    def select_channel_remove_stripes(self, ch):
        """select one channel in live mode, with remove_stripes functionality"""
        if ch == 488:
            self._select_channel_remove_stripes(self.ch_dio0)
        elif ch == 561:
            self._select_channel_remove_stripes(self.ch_dio1)
        else:
            raise ValueError('Channel not supported!')

    def _set_up_retriggerable_counter(self, counter):
        """set up a retriggerable counter task"""
        task_ctr = nidaqmx.Task()
        task_ctr.co_channels.add_co_pulse_chan_freq(counter,
                                                    idle_state=nidaqmx.constants.Level.LOW,
                                                    freq=self.sampling_rate)
        task_ctr.timing.cfg_implicit_timing(sample_mode=nidaqmx.constants.AcquisitionType.FINITE,
                                            samps_per_chan=self.num_samples)
        task_ctr.triggers.start_trigger.cfg_dig_edge_start_trig(trigger_source=self.PFI0,
                                                                trigger_edge=nidaqmx.constants.Slope.RISING)
        task_ctr.triggers.start_trigger.retriggerable = True
        return task_ctr

    def _retriggable_task(self, task, data):
        """set up a retriggeable task using a counter
        Using counter0 as default"""
        task.timing.cfg_samp_clk_timing(rate=self.sampling_rate,
                                        source=self.ch_ctr0_internal_output,
                                        sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)

        task_ctr = self._set_up_retriggerable_counter(self.ch_ctr0)

        task.write(data)
        task.start()
        task_ctr.start()
        print("type s to abort:")
        while input() != "s":
            time.sleep(0.05)

        task.stop()
        task_ctr.stop()
        task_ctr.close()
        task.close()

    def _retriggable_task_remove_stripes(self, task_do, data_do, task_ao, data_ao):
        """set up a retriggeable task using a counter
        Using counter0 as default"""
        task_do.timing.cfg_samp_clk_timing(rate=self.sampling_rate,
                                           source=self.ch_ctr0_internal_output,
                                           sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)
        task_ao.timing.cfg_samp_clk_timing(rate=self.sampling_rate,
                                           source=self.ch_ctr0_internal_output,
                                           sample_mode=nidaqmx.constants.AcquisitionType.CONTINUOUS)

        task_ctr = self._set_up_retriggerable_counter(self.ch_ctr0)

        task_do.write(data_do)
        task_ao.write(data_ao)
        task_do.start()
        task_ao.start()
        task_ctr.start()
        print("type s to abort:")
        while input() != "s":
            time.sleep(0.05)

        task_do.stop()
        task_ao.stop()
        task_ctr.stop()
        task_ctr.close()
        task_do.close()
        task_ao.close()


if __name__ == "__main__":
    daq_card = NIDaq(
            exposure=0.020,
            nb_timepoints=600,
            # scan_step=4,
            scan_step=0.310 * 2,  # TTL100
            # scan_step=0.376 * 4,  # TTL165
            # scan_step=0.103 * 2,    # TTL165
            stage_scan_range=2500,
            vertical_pixels=1024,  # used to calculate the readout time
            offset_view1=1580,
            offset_view2=1720,
            # offset_view1=1650,
            # offset_view2=1650,
            # view1_galvo1=4.52,
            # view1_galvo2=-4.05,
            # view2_galvo1=-4.18,
            # view2_galvo2=3.90,
            view1_galvo1=4.52,
            view1_galvo2=-4.00,
            view2_galvo1=-4.18,
            view2_galvo2=4.00,
            stripe_reduction_range=0.3,
            stripe_reduction_offset=-0.58
    )
    """Methods are separated into live mode and aacquisition mode"""

    # time lapse mode

    # live mode
    # daq_card.select_view(2)
    # daq_card.select_channel_remove_stripes(488)

    # daq_card.select_channel(488)

    # for 561
    daq_card.select_view(1)
    daq_card.select_channel_remove_stripes(488)

    # daq_card.select_channel(561)

    # daq_card.acquire_stacks(channels=[488], view=0)
    # daq_card.acquire_stacks(channels=[561], view=2)
    # daq_card.acquire_stacks(channels=[488, 561], view=0)

