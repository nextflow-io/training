// Courses shown on the homepage, in the course picker and in My courses.
// Titles, parts and durations come from the pages themselves (`duration:` frontmatter on each
// numbered lesson); this file only holds what the pages can't say.
//
// To add a course: add its pages to sidebar.mjs, then add an entry here.
// `id` is the course folder (its overview page's slug). The order here is the order learners see.
// `icon` is a name from src/icons.mjs. `buildsOn` lists the courses to take first (empty: a good first course).

export const groups = {
  run: 'Run pipelines',
  write: 'Write pipelines',
  science: 'Nextflow for Science',
};

export const courses = [
  {
    id: 'nextflow_run',
    buildsOn: [],
    group: 'run',
    icon: 'play',
    blurb: 'Launch pipelines, find outputs, use containers.',
    level: 'Beginner',
    next: ['nfcore_use', 'config_exec'],
  },
  {
    id: 'hello_nextflow',
    buildsOn: [],
    group: 'write',
    icon: 'code',
    blurb: 'Processes, channels, modules, containers and config.',
    level: 'Beginner',
    next: ['nfcore_build', 'side_quests'],
  },
  {
    id: 'config_exec',
    buildsOn: ['nextflow_run'],
    group: 'run',
    icon: 'settings',
    blurb: 'Resources, retries and execution profiles.',
    level: 'Intermediate',
    next: ['nfcore_use', 'seqera_scale'],
  },
  {
    id: 'nfcore_use',
    buildsOn: ['nextflow_run'],
    group: 'run',
    icon: 'package',
    blurb: 'Find, run and configure community pipelines.',
    level: 'Beginner',
    next: ['seqera_scale', 'nfcore_build'],
  },
  {
    id: 'nfcore_build',
    buildsOn: ['hello_nextflow', 'nfcore_use'],
    group: 'write',
    icon: 'hammer',
    blurb: 'Template, modules and input validation.',
    level: 'Intermediate',
    next: ['side_quests'],
  },
  {
    id: 'seqera_scale',
    buildsOn: ['nextflow_run', 'nfcore_use'],
    group: 'run',
    icon: 'cloud',
    blurb: 'Launch and monitor runs in Seqera Platform.',
    level: 'Intermediate',
    next: ['config_exec'],
  },
  {
    id: 'nf4_science/genomics',
    buildsOn: ['hello_nextflow'],
    group: 'science',
    icon: 'dna',
    blurb: 'Develop a variant calling pipeline, from per-sample calls to joint calling on a cohort.',
    level: 'Intermediate',
    next: ['side_quests', 'nfcore_build'],
  },
  {
    id: 'nf4_science/rnaseq',
    buildsOn: ['hello_nextflow'],
    group: 'science',
    icon: 'chart',
    blurb: 'Develop a bulk RNAseq processing pipeline, from a single sample to multi-sample paired-end data.',
    level: 'Intermediate',
    next: ['side_quests', 'nfcore_build'],
  },
  {
    id: 'nf4_science/imaging',
    buildsOn: [],
    group: 'science',
    icon: 'image',
    blurb: 'Run and configure imaging pipelines, using nf-core/molkart as the worked example.',
    level: 'Intermediate',
    next: ['nfcore_use', 'config_exec'],
  },
];

// Side Quests are listed from the "Side Quests" group in sidebar.mjs; this is its overview page.
export const sideQuests = { id: 'side_quests', group: 'write', icon: 'compass' };

// The homepage's goal tabs: which courses each one shows, in order.
export const goals = {
  run: ['nextflow_run', 'config_exec', 'nfcore_use', 'seqera_scale'],
  write: ['hello_nextflow', 'nfcore_use', 'nfcore_build'],
};

// Shown on course overview pages. Keep in step with courseInfo.defaultTechnicalRequirements.
export const nextflowVersion = '25.10.2';
