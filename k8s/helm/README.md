# rcsb-chem-search

K8s rcsb-chem-search related yaml files

## Setting up persistent volumes

The rcsb-chem-search persistent volume is used by the client application (defined in `filesystem.yaml`);
it is also updated by the ETL loader pods during the weekly workflow.
Independent persistent volumes and PVC are present for individual paths explained below.

### A-B Pathing Strategy

This service is deployed using the A-B pathing strategy for maintaining the weekly update cadence.
See the https://github.com/rcsb/devops-k8s-path-operator repository for additional details.
