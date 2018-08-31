# -*- coding: utf-8 -*-
"""
Created on Tue Mar  6 10:16:55 2018

@author: jdo

SpecsLab Prodigy uses "sle" file format, which is a SQLite database.
use http://inloop.github.io/sqlite-viewer/ to browse a file.




File Layout
^^^^^^^^^^^

the metadata we need is in the 'Schedule' XML String of the 'Configuration'
table

SELECT * FROM 'Configuration' LIMIT 0,30

Key             Value
Version	        1.11
AppVersion      4.37.2-r68309
Schedule        <!DOCTYPE SpecsLabSchedule><Experiment Version="1.1">...
Context         <Context><DeviceContext><DeviceMapping>...
MetaInfo        <MetaInfo/>

the regions have IDs, which can be found in the <Spectrum > tag of the Schedule
tree

<Spectrum ID="4" Name="Sb 3d"> <Comment></Comment> ... </Spectrum>

the count rates are in the table CountRateData

SELECT * FROM 'CountRateData' WHERE RawID=4

RawID   Offset  NextOffset  ChunkSize   Data
4       0       1218        30450       ...bytes (8 byte is one double)

for an MCD detector, the count rate is interlaced, so every (N+i*25)th value
belongs to channel N

"""

import sqlite3
import array
import re
import datetime
from xml.etree import ElementTree as ET
import logging
import numpy as np
from XPS.RegionBase import RegionBase
from XPS.cached_property import cached_property

_LOG = logging.getLogger(__name__)
# Make the logger follow the logging setup from the caller
_LOG.addHandler(logging.NullHandler())


class SpecsSLEFile(list):
    """This is the top structure for a parsed file which represents a list
    of RegionGroups

    The class contains a 'filepath' attribute.

    """

    def __init__(self, filepath, encoding=None):
        """Parse the XML and initialize the internal variables"""
        super(SpecsSLEFile, self).__init__()

        conn = sqlite3.connect(filepath)
        cursor = conn.cursor()

        cursor.execute("SELECT * FROM 'Configuration' WHERE Key='Schedule'")
        temp = cursor.fetchone()
        xmlstring = temp[1]

        if encoding:
            root = ET.fromstring(xmlstring.encode('utf-8'))
        else:
            try:
                root = ET.fromstring("<root>" + xmlstring[53:-13] + "</root>")
            except ET.ParseError:
                print('#####\nParsing of the XML file failed. Possibly the '
                      'XML is mal-formed or you need to supply the encoding '
                      'of the XML file.\n\n###Traceback:')
                raise

        for element in root.findall('Schedule/XPS/SpectrumGroup'):
            _LOG.debug('Found region group: %s', element)
            self.append(RegionGroup(element, cursor))

        del(cursor)
        conn.close()

    def __repr__(self):
        """Returns class representation"""
        return '<{}>'.format(self.__class__.__name__)

    def __str__(self):
        """Returns str representation"""
        out = self.__repr__()
        for region_group in self:
            for line in region_group.__str__().split('\n'):
                out += '\n    ' + line
        return out




class RegionGroup(list):
    """Class that represents a region group, which consist of a list of
    regions

    The class contains a 'name' and and 'parameters' attribute.

    """

    def __init__(self, xml, cursor):
        """Initializes the region group

        Expects to find 3 subelement; the name, regions and
        parameters. Anything else raises an exception.

        Parsing parameters is not supported and therefore logs a
        warning if there are any.

        """
        super(RegionGroup, self).__init__()

        # Get name, find a string tag with attribute 'name' with value 'name'
        self.name = xml.attrib['Name']
        self.node_id = xml.attrib['ID']

        # the RegionGroup contains the MCD detector calibration
        cursor = cursor
        self.detectors = []
        cursor.execute("SELECT * FROM 'NodeData' WHERE Node={}".format(self.node_id))
        temp = cursor.fetchone()
        elements = ET.fromstring(temp[1]).findall('./DetectorCalibration/Detector')
        for element in elements:
            gain = float(element.attrib['Gain'])
            position = float(element.attrib['Position'])
            shift = float(element.attrib['Shift'])
            self.detectors.append({'gain': gain, 'position': position, 'shift': shift})

        # and it knows the Spectra it contains
        for element in xml.findall('Spectrum'):
            _LOG.debug('Found region: %s', element)
            self.append(Region(element, cursor, self))

    def __repr__(self):
        """Returns class representation"""
        return '<{}(name=\'{}\')>'.format(self.__class__.__name__, self.name)

    def __str__(self):
        """Return the class str representation"""
        out = self.__repr__()
        for region in self:
            out += '\n    ' + region.__str__()
        return out


class Region(RegionBase):
    """Class that represents a region

    The class contains attributes for the items listed in the
    'information_names' class variable.

    Some useful ones are:
    * **name**: The name of the region
    * **ID**: The ID it goes by in the SQLite file
    * **region**: Contains information like, DwellTime, AnalysisMethod,
        ScanDelta, ExcitationEnergy etc.
        Note that these are slightly different than in the SpecsLab 2 XML
        format

    All auxiliary information is also available from the 'info'
    attribute.

    """
    # pylint: disable=too-many-instance-attributes
    # 17 is reasonable in this case.

    information_names = ['AnalysisMethod', 'AnalyzerLensMode',
                         'AnalyzerLensVoltage', 'AnalyzerSlit', 'BiasVoltage',
                         'CurvesPerScan', 'DetectorVoltage', 'DwellTime',
                         'ExcitationEnergy', 'KineticEnergy',
                         'KineticEnergyBase', 'PassEnergy', 'ScanDelta',
                         'ScanMode', 'ValuesPerCurve', 'Workfunction']

    def __init__(self, xml, cursor, parent):
        """Parse the XML and initialize internal variables
        Args:
            xml (xml.etree.ElementTree.Element): The region XML element
        """
        # Parse information items
        self.region = {}
        self.name = xml.attrib['Name']
        self.node_id = xml.attrib['ID']
        self.parent = parent

        cursor.execute("SELECT * FROM 'NodeData' WHERE Node={}".format(self.node_id))
        temp = cursor.fetchone()
        element = ET.fromstring(temp[1]).find('./AnalyzerSpectrumParameters')
        if element is not None:
            for name in self.information_names:
                nom = self.convert(name)
                self.region[nom] = element.attrib[name]
                try:
                    self.region[nom] = float(element.attrib[name])
                except ValueError:
                    pass
                # Dynamically create attributes for all the items
                setattr(self, name, self.region[nom])

        cursor.execute("SELECT * FROM 'Spectrum' WHERE Node={}".format(self.node_id))
        spec = cursor.fetchone()
        if spec is not None:
            self.mcd_head = spec[4]
            self.mcd_tail = spec[5]
            self.mcd_num = spec[1]
            self.num_samples = spec[3]
            self.dim_transm = spec[7]
            self.num_transm = spec[8]
            self.transm = array.array('d', spec[9])
            
            # be compatible with the XML version:
            self.info = {'transmission': np.reshape(np.array(self.transm),
                         (self.dim_transm, self.num_transm))}

            self.detectors = self.parent.detectors

        # get the RawID where the counts are stored that belong to this NodeID
        cursor.execute("SELECT * FROM 'RawData' WHERE Node={}".format(self.node_id))
        temp = cursor.fetchone()
        if temp is not None:
            self.raw_id = temp[0]
            self.scan_date = temp[2]
            self.unix_timestamp = datetime.datetime.strptime(temp[2],
                                                             "%Y-%b-%d %H:%M:%S.%f"
                                                            ).timestamp()

            # now use the RawID to get the counts
            cursor.execute("SELECT * FROM 'CountRateData' WHERE RawID={}".format(self.raw_id))
            temp = cursor.fetchone()
            if temp is not None:
                self.y_avg_counts_data = np.array(array.array('d', temp[4]))

            _LOG.debug('Creating y_avg_counts values from scans')

    @staticmethod
    def convert(name):
        """converts the CamelCase names in the SQLite file to cunder_score_format"""
        str_1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
        return re.sub('([a-z0-9])([A-Z])', r'\1_\2', str_1).lower()

    @cached_property
    def y_avg_counts(self):
        """Returns the average counts from a cycle as a Numpy array"""
        return self.y_avg_counts_data
