from examples.A0_SimpleExample import main as run_simple_example
from examples.A1_SimpleReactionExample import main as run_reaction_example
from examples.A2_SimpleReactionExample_StepByStepSim import main as run_step_by_step_example
from examples.B0_SimplePackedBedExample import main as run_packed_bed_example
from examples.B1_SimplePackedBedExample_StepByStepSim import main as run_packed_bed_step_example
from examples.B2_PackedBedExample_StoppingFlowRate import main as run_packed_bed_stopping_flow_example
from examples.C1_LoadAndPrepareMeasurements import main as load_and_prepare_measurements

if __name__ == "__main__":
    #run_simple_example()
    #run_reaction_example()
    #run_step_by_step_example()

    #run_packed_bed_example()
    #run_packed_bed_step_example()
    #run_packed_bed_stopping_flow_example()

    load_and_prepare_measurements() 

    input("Press Enter to exit...")