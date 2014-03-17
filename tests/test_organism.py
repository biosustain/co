import unittest
from colib import Component
from colib.mutations import SNP, Mutation
from colib.organisms import HaploidOrganism
from colib.position import Range
from colib.storage import Storage


class HaploidOrganismTestCase(unittest.TestCase):

    def setUp(self):
        self.storage = storage = Storage()
        storage.components.add(Component('genome-1'))
        storage.components.add(Component('part-1'))
        storage.components.add(Component('plasmid-1'))

        strain = HaploidOrganism('strain-1')
        # strain.contigs.add(Contig('genome', self.storage.components['genome-1']))
        # strain.contigs.add(Contig('plasmid-a', self.storage.components['plasmid-1']))

        storage.add_organism(strain)

    def tearDown(self):
        pass

    def test_create_strain(self):
        strain = self.storage.organisms['strain-1']
        self.assertEqual(strain.as_dict(), {'display_id': 'strain-1', 'contigs': [{'name': 'genome', 'component': 'genome-1'}]})

    def test_mutate_strain(self):
        strain = self.storage.organisms['strain-1']
        feature_4 = strain.contigs[0].features[4]
        feature_5 = strain.contigs[0].features[5]

        mutations = [
            SNP(100, 'A'),
            SNP(105, 'G'),
            Mutation(Range(40, 45), None),
            Mutation(Range(20, 21), 'GATGA'),
            Mutation(Range(feature_4.start, feature_4.end + 20), None),
            Mutation(Range(feature_5.start, feature_5.start), self.storage.components['part-1'])
        ]

        new_genome = strain.contigs[0].mutate(mutations)

        new_strain = HaploidOrganism('strain-1-1', parents=[strain])
        new_strain.contigs.add(new_genome) # replaces genome.

        self.storage.add_organism(new_strain) # saves strain.

        self.assertEqual(False, feature_4 in new_strain)
        self.assertEqual(True, feature_5 in new_strain)

        self.assertEqual(len(new_strain.diff(strain)), 2)


    def test_list_features(self):
        "Expect all features from genome-1 and plasmid-1."
        pass