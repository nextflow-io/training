import { defineCollection } from 'astro:content';
import { z } from 'astro/zod';
import { docsLoader, i18nLoader } from '@astrojs/starlight/loaders';
import { docsSchema, i18nSchema } from '@astrojs/starlight/schema';

const additionalInformation = z
  .object({
    technical_requirements: z.union([z.literal(true), z.string()]).optional(),
    learning_objectives: z.array(z.string()).optional(),
    audience_prerequisites: z.array(z.string()).optional(),
    videos: z.string().optional(),
    videos_playlist: z.string().url().optional(),
  })
  .refine((v) => !(v.videos && v.videos_playlist), 'videos and videos_playlist are mutually exclusive');

export const collections = {
  docs: defineCollection({
    loader: docsLoader(),
    schema: docsSchema({
      extend: z.object({
        page_type: z.literal('index_page').optional(),
        index_type: z.enum(['course', 'module']).optional(),
        additional_information: additionalInformation.optional(),
      }),
    }),
  }),
  i18n: defineCollection({
    loader: i18nLoader(),
    schema: i18nSchema({
      extend: z
        .object({
          'courseInfo.courseSummary': z.string(),
          'courseInfo.additionalInformation': z.string(),
          'courseInfo.technicalRequirements': z.string(),
          'courseInfo.learningObjectives': z.string(),
          'courseInfo.audiencePrerequisites': z.string(),
          'courseInfo.courseVideos': z.string(),
          'courseInfo.defaultTechnicalRequirements': z.string(),
          'courseInfo.defaultVideos': z.string(),
          'courseInfo.viewPlaylist': z.string(),
          'consent.title': z.string(),
          'consent.description': z.string(),
          'consent.posthog': z.string(),
        })
        .partial(),
    }),
  }),
};
