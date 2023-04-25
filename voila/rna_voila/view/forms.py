import numpy as np

from flask_wtf import FlaskForm
from wtforms import BooleanField, StringField, SelectField
from wtforms.fields.html5 import DecimalRangeField

class LsvFiltersForm(FlaskForm):
    a5ss = BooleanField('5 Prime')
    a3ss = BooleanField('3 Prime')
    exon_skipping = BooleanField('Exon Skipping')
    source = BooleanField('Source')
    target = BooleanField('Target')
    binary = BooleanField('Binary')
    complex = BooleanField('Complex')
    intron_retention = BooleanField("Intron Retention")


class DeltaPsiFiltersForm(FlaskForm):
    dpsi_threshold = StringField('abs(E(dPSI)) Threshold', default=0.2)
    confidence_threshold = DecimalRangeField('Confidence Threshold', default=0.9)

class HeterogenFiltersForm(FlaskForm):
    dpsi_threshold = StringField('abs(E(dPSI)) Threshold', default=0.2)
    stat_threshold = DecimalRangeField('P-Value Threshold', default=0.05)
    stat_type = SelectField('P-Value stat', coerce=int, choices=[])

    def __init__(self, choices):
        super().__init__()

        for i, choice in enumerate(choices):
            self.stat_type.choices.append((i, choice))
