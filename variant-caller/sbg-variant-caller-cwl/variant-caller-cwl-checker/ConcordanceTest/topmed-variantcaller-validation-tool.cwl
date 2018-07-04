{
    "class": "CommandLineTool",
    "cwlVersion": "v1.0",
    "id": "vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/10",
    "baseCommand": [],
    "inputs": [
        {
            "sbg:category": "Input Files",
            "id": "Expected_result",
            "type": "File",
            "label": "Expected result",
            "sbg:fileTypes": "TAR.GZ"
        },
        {
            "sbg:category": "Input Files",
            "id": "Workflow_result",
            "type": "File",
            "sbg:fileTypes": "TAR.GZ"
        }
    ],
    "outputs": [
        {
            "id": "Match",
            "type": "string",
            "outputBinding": {
                "loadContents": true,
                "glob": "match.vcf",
                "outputEval": "${\n    if (self[0].size) {\n        return \"Failure\"\n    } else {\n        return \"Success\"\n    }\n}"
            }
        }
    ],
    "label": "Validation_tool",
    "arguments": [
        {
            "position": 0,
            "prefix": "",
            "shellQuote": false,
            "valueFrom": "${\n    var comm =\"\"\n    comm += \"mkdir workflow_result expected_result\"\n    comm += \" && tar -xzvf \" + inputs.Workflow_result.path + \" --strip-components 1 -C workflow_result/\"\n    comm += \" && tar -xzvf \" + inputs.Expected_result.path + \" --strip-components 1 -C expected_result/\"\n    comm += \" && mkdir match\"\n    comm += \" && bcftools isec -p match expected_result/*.vcf.gz workflow_result/*.vcf.gz \"\n    comm += \" && vcf-concat match/0000.vcf match/0001.vcf >> match2.vcf\"\n    comm += \" && sed '/^#/d' match2.vcf >> match.vcf\"\n    \n    return comm\n    \n}"
        }
    ],
    "requirements": [
        {
            "class": "ShellCommandRequirement"
        },
        {
            "class": "ResourceRequirement",
            "ramMin": 100,
            "coresMin": 1
        },
        {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/vladimir_obucina/topmed:validation_tool"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "sbg:modifiedOn": 1527588853,
    "sbg:latestRevision": 10,
    "sbg:sbgMaintained": false,
    "sbg:validationErrors": [],
    "sbg:id": "vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/10",
    "sbg:createdOn": 1527513179,
    "sbg:contributors": [
        "vladimir_obucina"
    ],
    "sbg:image_url": "https://igor.sbgenomics.com/ns/brood/images/vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline/validation_tool/10.png",
    "sbg:revisionsInfo": [
        {
            "sbg:revision": 0,
            "sbg:modifiedOn": 1527513179,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": null
        },
        {
            "sbg:revision": 1,
            "sbg:modifiedOn": 1527522399,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": "First version with bcftools and vcftools"
        },
        {
            "sbg:revision": 2,
            "sbg:modifiedOn": 1527522451,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 3,
            "sbg:modifiedOn": 1527522628,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 4,
            "sbg:modifiedOn": 1527522900,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 5,
            "sbg:modifiedOn": 1527523085,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 6,
            "sbg:modifiedOn": 1527585970,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 7,
            "sbg:modifiedOn": 1527587057,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": "UPDATE: working locally"
        },
        {
            "sbg:revision": 8,
            "sbg:modifiedOn": 1527587374,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 9,
            "sbg:modifiedOn": 1527588145,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        },
        {
            "sbg:revision": 10,
            "sbg:modifiedOn": 1527588853,
            "sbg:modifiedBy": "vladimir_obucina",
            "sbg:revisionNotes": ""
        }
    ],
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "sbg:modifiedBy": "vladimir_obucina",
    "sbg:publisher": "sbg",
    "sbg:projectName": "TOPMed Freeze 3a Variant Calling Pipeline",
    "sbg:revisionNotes": "",
    "sbg:appVersion": [
        "v1.0"
    ],
    "sbg:createdBy": "vladimir_obucina",
    "sbg:project": "vladimir_obucina/topmed-freeze-3a-variant-calling-pipeline",
    "sbg:revision": 10
}