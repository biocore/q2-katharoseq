import qiime2.plugin.model as model


STATS_HEADER = ['sample-id',
                'log_total_reads',
                'estimated_biomass_per_pcrrxn',
                'estimated_biomass_per_dnarxn',
                'extraction_mass_g',
                'estimated_cells_per_g',
                'log_estimated_cells_per_g']


class EstimatedBiomassFmt(model.TextFileFormat):
    def sniff(self):
        line = open(str(self)).readline()
        hdr = line.strip().split(',')

        return hdr == STATS_HEADER

    def validate(self, *args):
        pass


EstimatedBiomassDirFmt = model.SingleFileDirectoryFormat(
    'EstimatedBiomassDirFmt', 'est_biomass.csv', EstimatedBiomassFmt)
