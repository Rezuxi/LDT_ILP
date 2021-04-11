# LDT_ILP

ILP solver for LDT editing with the least amount of edits.

The results along with the unedited graph and runtime is stored in json format.

The filenames of the exact results include the number of nodes, n (int), if it's only adding/deleting edges (0 or 1) and an ID (int).

filename + " _ n _ (only_add) _ (only_delete) _ ID.json"