import abc

# ###############################
# Data loading and Boilerplate #
################################


def pipeline_run(pipeline):
    """ Exception catching for all the pipelines """
    try:
        return pipeline.run()
    except KeyboardInterrupt:
        if pipeline.logger:
            pipeline.logger.info("MAJIQ manually interrupted. Avada kedavra...")


class BasicPipeline:
    def __init__(self, args):
        """Basic configuration shared by all pipelines"""

        self.__dict__.update(args.__dict__)
        self.nthreads = args.nproc
        self.dpsi = False

    @abc.abstractmethod
    def run(self):
        """This is the entry point for all pipelines"""
        return
