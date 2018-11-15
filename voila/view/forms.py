import numpy as np

from flask_wtf import FlaskForm
from wtforms import BooleanField, StringField, SelectField


class LsvFiltersForm(FlaskForm):
    a5ss = BooleanField('5 Prime')
    a3ss = BooleanField('3 Prime')
    exon_skipping = BooleanField('Exon Skipping')
    source = BooleanField('Source')
    target = BooleanField('Target')
    binary = BooleanField('Binary')
    complex = BooleanField('Complex')


class DeltaPsiFiltersForm(FlaskForm):
    dpsi_threshold = StringField('|E(dPSI)| Threshold', default=0.2)
    confidence_threshold = StringField('Confidence Threshold', default=0.95)
    # confidence_threshold = SelectField('Confidence Threshold', choices=list((x, x) for x in np.linspace(0, 1, 10)))
