## Agent skills

### Issue tracker

Issues are tracked in this repository's GitHub Issues. See `docs/agents/issue-tracker.md`.

### Triage labels

The default five canonical triage labels are used. See `docs/agents/triage-labels.md`.

### Domain docs

This repository uses a single-context domain-doc layout. See `docs/agents/domain.md`.

# Remote Compute and Slurm Execution Standard

This section defines the mandatory rules for interacting with remote compute servers and Slurm clusters.

These rules are project-independent.

Project-specific values such as SSH hosts, remote paths, environments, partitions, accounts, datasets, and resource requirements must be defined in project configuration or project scripts, not hard-coded into this policy.

------

## 1. Core Architecture

Assume the following general execution model:

```text
Local development machine
        |
        | SSH / rsync
        v
Optional SSH jump host
        |
        v
Remote login node
        |
        | Slurm
        v
Compute nodes
```

The local machine is the development environment.

The remote cluster is a compute environment.

An SSH jump host, VPN gateway, or ProxyJump configuration is infrastructure and should be treated as transparent.

Use the SSH host aliases already defined by the user.

Do not modify:

- VPN settings;
- network routes;
- SSH jump-host topology;
- firewall configuration;
- system-wide SSH configuration;

unless explicitly asked to do so.

------

# 2. Source-of-Truth Rule

The local repository is the authoritative source of code.

Always:

- inspect code locally;
- edit code locally;
- review diffs locally;
- perform Git operations locally unless explicitly required otherwise.

Do not edit project source code directly on:

- the SSH jump host;
- the remote login node;
- compute nodes;
- temporary Slurm working directories.

Remote source trees are execution copies or immutable run snapshots.

Any code fix discovered from remote execution must be applied to the local repository first and then redeployed.

------

# 3. Never Hard-Code Remote Infrastructure

Do not hard-code project-specific infrastructure into source code or this policy.

Values such as:

```text
SSH host
remote deployment directory
remote run directory
Conda or virtual environment
Slurm account
Slurm partition
CPU count
memory
wall time
GPU requirements
```

must come from project configuration, environment files, wrapper scripts, or Slurm scripts.

Prefer a project-local configuration such as:

```text
dev/remote.env
```

or an equivalent project-defined configuration.

Do not guess missing infrastructure values.

------

# 4. Preferred Remote Interface

When the project provides standardized remote helper commands, use them instead of constructing ad-hoc SSH, rsync, or Slurm commands.

Preferred interfaces are:

```bash
./dev/push
./dev/remote-test
./dev/slurm-submit
./dev/slurm-status
./dev/slurm-log
```

Their conceptual responsibilities are:

```text
dev/push
    synchronize the current local development tree

dev/remote-test
    run short remote validation

dev/slurm-submit
    submit reproducible compute jobs

dev/slurm-status
    inspect scheduler/accounting state

dev/slurm-log
    inspect stdout/stderr
```

If equivalent project-specific wrappers exist under different names, prefer those.

Do not bypass a project-provided wrapper without a clear reason.

------

# 5. Local vs Remote Decision Rule

Before executing a command, classify it.

## Run locally when the task is primarily:

- source inspection;
- editing;
- formatting;
- linting;
- type checking;
- Git inspection;
- small unit testing;
- fixture-based testing;
- small synthetic-data testing;
- lightweight preprocessing;
- configuration validation.

## Use the remote environment when the task requires:

- software available only on the cluster;
- remote datasets;
- cluster-specific environments;
- integration with remote infrastructure;
- short environment-dependent smoke tests.

## Use Slurm when the task is:

- computationally expensive;
- long-running;
- CPU-intensive;
- memory-intensive;
- GPU-dependent;
- parallel;
- distributed;
- data-intensive;
- a large simulation;
- a parameter sweep;
- a replicate campaign;
- a full validation run;
- a large model-training or inference task.

When uncertain whether substantial computation should run on the login node or under Slurm:

**use Slurm.**

------

# 6. Login Node Policy

Treat the remote login node as an orchestration environment, not a compute node.

Appropriate login-node operations include:

- file inspection;
- environment inspection;
- checking installed software;
- lightweight metadata queries;
- short smoke tests;
- small diagnostic commands;
- `rsync`;
- `git` inspection when required;
- `sbatch`;
- `squeue`;
- `sacct`;
- log inspection;
- result metadata inspection.

Do not run substantial scientific or computational workloads directly on the login node.

Do not launch long-running workloads using:

```bash
nohup ...
command &
screen ...
tmux ...
```

as substitutes for Slurm scheduling.

Interactive compute work should use the cluster's supported Slurm interactive mechanism, such as `srun` or `salloc`, when appropriate and when supported by the cluster.

------

# 7. Synchronization Rule

Before remote testing, ensure the remote development copy corresponds to the intended local source state.

Prefer:

```bash
./dev/push
```

or the project's equivalent synchronization command.

Synchronization must be unidirectional:

```text
local source
    ->
remote deployment copy
```

Do not use bidirectional synchronization for authoritative source code.

Do not treat shared cloud folders, network drives, or automatically synchronized directories as proof that the remote execution copy is current.

Explicit deployment is preferred.

------

# 8. Protect Remote Data

Source deployment and scientific result storage must be separate concerns.

A disposable deployment directory may be overwritten during development.

Result directories, datasets, run archives, and provenance directories must not be treated as disposable.

Never run destructive commands such as:

```bash
rm -rf
rsync --delete
find ... -delete
```

against remote data, result, run, or dataset directories unless the target has been explicitly verified as disposable.

Never delete previous experimental results merely to make a new run succeed.

If the destination is ambiguous, do not perform destructive synchronization.

------

# 9. Remote Testing Policy

Remote tests should validate the environment before expensive computation.

Preferred sequence:

```text
local modification
    ->
local lightweight test
    ->
remote synchronization
    ->
remote smoke/integration test
    ->
Slurm submission
```

Use remote tests for short validation only.

Examples include:

- imports;
- environment checks;
- small integration tests;
- one small input;
- one small replicate;
- configuration loading;
- dependency checks.

Do not turn `remote-test` into a substitute for Slurm.

------

# 10. Slurm Is the Execution Boundary for Heavy Work

Heavy cluster computation must be submitted through Slurm.

Prefer:

```bash
./dev/slurm-submit <job-script>
```

when available.

Otherwise use the project's approved Slurm submission mechanism.

Avoid directly running substantial workloads using:

```bash
ssh remote 'python large_job.py'
```

or equivalent commands.

The SSH connection should orchestrate the job, not remain attached to the entire computation.

A successful submission is not the same as a successful computation.

------

# 11. Immutable Run Principle

A submitted Slurm job must execute a stable source version.

Do not allow queued or running jobs to depend on a mutable development directory that may later be overwritten by another deployment.

Preferred architecture:

```text
remote-deploy/
    mutable development copy

remote-runs/
    RUN_ID/
        source/
        provenance/
        results/
        logs/
```

For substantive runs, create an immutable or effectively immutable source snapshot at submission time.

Subsequent local edits or `dev/push` operations must not change the source tree used by already-submitted jobs.

------

# 12. Provenance Requirement

Every substantive Slurm submission should preserve enough information to identify exactly what was submitted.

Where available, record:

```text
run ID
Slurm job ID
submission timestamp
Git commit
Git branch
Git dirty status
Git diff when dirty
source content hash
Slurm script
Slurm script arguments
relevant execution configuration
remote run directory
```

For reproducibility-critical work, also preserve when practical:

```text
software environment
dependency versions
container image or environment identifier
random seed
replicate identifier
input manifest
parameter/configuration hash
resource request
```

The run must remain traceable even if the local repository changes later.

------

# 13. Git Dirty-State Policy

A dirty Git working tree does not automatically prohibit exploratory computation.

For exploratory or development work:

```text
dirty working tree
    ->
allowed only if the exact dirty state is captured
```

Capture at minimum:

```text
Git commit
dirty=true
git diff
source snapshot or content hash
```

For formal, publication-quality, benchmark, production, or otherwise reproducibility-critical runs:

prefer a clean committed working tree.

If the project configuration requires a clean tree, respect that requirement and do not bypass it.

Never report a dirty run as though it corresponded exactly to its Git commit.

------

# 14. Slurm Submission Checks

Before submitting a substantial job:

1. confirm the intended source state;
2. run the smallest appropriate validation;
3. ensure required files are available;
4. inspect the Slurm script;
5. inspect important resource requests;
6. verify output/run paths;
7. submit through the standard wrapper;
8. capture the returned job ID.

Important resource parameters may include:

```text
partition
account
nodes
tasks
CPUs
memory
GPU resources
wall time
array size
```

Do not silently change scientific or computational parameters simply to make a job easier to schedule.

------

# 15. Large Job and Job-Array Safety

Before submitting a potentially expensive job, large job array, or large replicate campaign, determine whether the user explicitly requested execution.

If the user explicitly asked to run or submit the computation, submission may proceed after validating the configuration.

If the user only asked to:

- modify code;
- prepare an experiment;
- create a Slurm script;
- inspect configuration;
- estimate resources;

do not automatically launch a large computation.

Do not accidentally convert a smoke test into hundreds or thousands of scheduled tasks.

Always inspect array bounds and replicate counts before submission.

------

# 16. Job Identification

After `sbatch`, always capture the Slurm job ID.

Never refer to an active experiment only by phrases such as:

```text
the latest job
the last run
the current simulation
```

when a job ID or run ID is available.

Prefer:

```text
RUN_ID
SLURM_JOB_ID
```

as stable identifiers.

For arrays, preserve both:

```text
array job ID
array task ID
```

when relevant.

------

# 17. Job Monitoring Policy

Do not interpret queue delay as job failure.

For active or pending jobs, use scheduler state tools such as:

```bash
squeue
```

For completed or historical jobs, use accounting information such as:

```bash
sacct
```

or project wrappers such as:

```bash
./dev/slurm-status <JOB_ID>
```

Inspect scheduler reasons before taking corrective action.

Examples of states or conditions that require different responses include:

```text
PENDING
RUNNING
COMPLETED
FAILED
CANCELLED
TIMEOUT
OUT_OF_MEMORY
NODE_FAIL
PREEMPTED
```

Do not automatically resubmit a job simply because it is pending.

------

# 18. No Duplicate Submission Without Diagnosis

Before resubmitting a failed, missing, or apparently stalled run:

1. inspect Slurm state;
2. inspect stdout;
3. inspect stderr;
4. inspect existing output files;
5. inspect run metadata;
6. determine whether an existing job is still active;
7. determine the failure mode.

Do not blindly submit duplicate jobs.

If a retry is necessary, preserve the relationship between the original run and retry when the project supports it.

------

# 19. Failure Diagnosis

For failed remote jobs, distinguish between:

```text
code failure
configuration failure
environment failure
missing input
filesystem/path failure
permission failure
scheduler failure
resource exhaustion
timeout
out-of-memory termination
node failure
scientific/model failure
```

Do not treat every Slurm failure as a code defect.

Apply code fixes locally.

Then:

```text
fix locally
    ->
local validation
    ->
redeploy
    ->
small remote validation
    ->
resubmit if justified
```

------

# 20. Resource Changes Must Be Explicit

Do not silently increase or decrease resources.

Changes to parameters such as:

```text
--time
--mem
--cpus-per-task
--nodes
--ntasks
--gres
--partition
--array
```

can materially change cost, scheduling behavior, runtime behavior, or scientific execution.

When correcting a resource-related failure, explain the proposed resource change.

Never reduce scientific workload parameters merely to turn a failed scientific run into a technically successful one unless the task explicitly calls for a smaller validation run.

------

# 21. Results Stay Remote by Default

Large raw outputs should remain on remote storage unless they are needed locally.

Prefer retrieving:

- compact summaries;
- metrics;
- logs;
- reports;
- small tables;
- selected plots;
- debugging artifacts.

Do not automatically transfer:

- very large result trees;
- large simulation outputs;
- large intermediate files;
- datasets already available remotely.

When local analysis requires remote outputs, retrieve only the necessary subset where practical.

------

# 22. Execution Success Is Not Scientific Success

Always distinguish:

```text
job submitted
job started
job completed
software exited successfully
validation passed
scientific conclusion supported
```

These are different states.

A Slurm job reaching `COMPLETED` only means that Slurm observed successful process completion.

It does not by itself establish:

- correctness;
- scientific validity;
- statistical significance;
- model adequacy;
- benchmark superiority;
- reproducibility.

Base scientific conclusions on actual outputs and predefined evaluation criteria.

------

# 23. Never Fabricate Remote State

Never claim that:

```text
a job was submitted
a job is running
a job completed
a test passed
a file exists remotely
a result was produced
```

unless that state was actually observed through the relevant command or output.

If remote access is unavailable, report that limitation instead of assuming remote state.

------

# 24. Long-Running Jobs

Do not keep an interactive Codex command or SSH connection open solely to wait for a long-running Slurm job.

For long-running work:

```text
submit
    ->
capture job ID
    ->
return control
```

Later status checks should use the job ID.

Do not poll excessively.

Use explicit status checks when needed.

------

# 25. Remote Command Safety

Before executing a remote command, consider whether it can:

- delete data;
- overwrite results;
- alter permissions;
- terminate jobs;
- cancel arrays;
- modify shared environments;
- modify shared datasets;
- consume substantial cluster resources;
- affect other users.

For potentially destructive or high-impact operations, verify the target and intent before execution.

Commands such as:

```bash
scancel
rm
rm -rf
chmod -R
chown -R
rsync --delete
```

require particular care.

Never use broad destructive patterns when a narrower operation is possible.

------

# 26. Shared Cluster Etiquette

Assume the cluster is shared infrastructure.

Do not intentionally bypass the scheduler for compute workloads.

Do not consume login-node resources with heavy computation.

Do not attempt to circumvent:

- scheduler policies;
- resource quotas;
- account restrictions;
- queue policies;
- access controls.

Respect the site's configured Slurm policies.

------

# 27. Standard Agent Workflow

For tasks involving remote computation, follow this default sequence:

```text
1. Inspect locally.
2. Modify locally.
3. Run the smallest useful local validation.
4. Inspect the local diff.
5. Deploy/synchronize if remote execution is needed.
6. Run a short remote smoke or integration test if appropriate.
7. Determine whether the real workload requires Slurm.
8. Create or verify an immutable run snapshot.
9. Record provenance.
10. Submit through Slurm.
11. Capture RUN_ID and JOB_ID.
12. Return control rather than waiting on long-running work.
13. Inspect scheduler state when requested or required.
14. Inspect logs and results after completion.
15. Diagnose failures before retrying.
16. Make fixes locally.
17. Repeat from the smallest useful validation level.
```

------

# 28. Decision Summary

Use this decision tree:

```text
Does the task involve editing source code?
    YES -> work locally.

Is the computation small and independent of the cluster?
    YES -> run locally.

Does it require remote software/data but finish quickly?
    YES -> short remote test may be appropriate.

Is it expensive, long-running, parallel, memory-heavy,
GPU-dependent, data-heavy, or scientifically substantive?
    YES -> use Slurm.

Is the login node being used for substantial computation?
    YES -> stop and move the workload to Slurm.

Is this a formal or reproducibility-critical experiment?
    YES -> use immutable source snapshot + provenance.

Has the job already been submitted?
    YES -> identify it using RUN_ID/JOB_ID and inspect its state;
           do not blindly resubmit.

Did a remote job reveal a code defect?
    YES -> fix locally, redeploy, validate, then resubmit.

Is a destructive or high-impact operation required?
    YES -> verify target and intent before executing.
```

------

# 29. Separation of Concerns

Keep these layers separate:

```text
AGENTS.md
    = behavioral policy

project remote configuration
    = infrastructure values

dev/push
    = deployment

dev/remote-test
    = short remote validation

dev/slurm-submit
    = reproducible job submission

dev/slurm-status
    = scheduler/accounting inspection

dev/slurm-log
    = log inspection

Slurm scripts
    = job-specific resources and commands

application code
    = scientific or computational logic
```

Do not move infrastructure-specific details into this policy.

Do not move scientific parameters into generic remote-execution tooling unless they are explicitly part of that project's configuration.

------

# 30. Final Principle

The default philosophy is:

```text
Develop locally.
Deploy explicitly.
Validate cheaply.
Compute through Slurm.
Snapshot important runs.
Record provenance.
Identify jobs explicitly.
Diagnose before retrying.
Protect remote data.
Keep infrastructure separate from project logic.
```
