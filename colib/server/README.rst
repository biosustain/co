Component-lib Server
====================

Standalone REST API and parts library based on Flask, Flask-SQLAlchemy, Flask-Presst and component-lib.


Commands
========

::

    GET /component
    GET /component/{id}

    {
        "parent": None,
        "name": None,
        "type": "plasmid",
        "uri": "/component/1"
        "meta": {}
    }

    GET /component/{id}/sequence
    Returns the nucleotide sequence of the component.

    GET /component/{id}/features

    [{
        "pos": 12,
        "size": 15,
        "name": None,
        "type": "gene",
        "qualifiers": {
            "db_xref": ["GI:6319249"],
        },
        "source": ["/component/12", False]
    }]

    GET /component/{id}/features/intersect ? start end include_inherited=True

    GET /component/{id}/features/find ? type=Any name=Any qualifiers={ ... } xref=Any

    POST /component/mutate {mutations}
    Creates a copy of the component with the mutations applied. Alternatively a component can be created using
    POST /component with "parent", "mutations" and "_deleted_features". This is suggested if any additional changes to the component's features
    are made. Internally, the server will replay the mutation algorithm to remove any features that should be removed.

    GET /component/{id}/diff {other: component} ? algorithm=False

    Returns a list of mutations that need to be applied to the component to get the
    sequence of the other component. If an algorithm is not given, this only works for components within the same
    lineage.

    GET /strain/{id}

    {
        "name": None,
        "parent": "/strain/1",
        "components": {
            "I": {"type": "chromosome", "component": "/component/1"},
            "II": {"type": "chromosome", "component": "/component/3"}
        },
        "uri": "/strain/2"
    }

    GET /strain/{id}/features/find ..
    GET /strain/{id}/features/count? type=any name=any qualifiers={ ... }

    GET /strain/{id}/diff {strain}
    Returns a list of added, removed, or changed components.

    {
        "added": [],
        "removed": [],
        "changed": ["II"]
    }

    GET /strain/{id}/lineage

    [
        "strain/2",
        "strain/1"
    ]

    GET /strain/tree {restrict: array of strains}
    [
        ["strain/1", [["strain/2", None], ["strain/3", None]],
        ["strain/5", [["strain/6", ["strain/7", None]]]
    ]


Lock
====

A component is locked for further mutation, or deletion, when any other components inherit from it. In future releases,
this lock may be overridden.

Component and feature types
===========================

All types should map to Sequence Ontology terms.