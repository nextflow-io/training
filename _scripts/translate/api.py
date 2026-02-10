"""Claude API interaction — async-first implementation."""

from __future__ import annotations

import asyncio
import random

import anthropic

from .config import (
    BASE_DELAY,
    MAX_CONTINUATIONS,
    MAX_RETRIES,
    MAX_TOKENS,
    MODEL,
    REQUEST_TIMEOUT,
    TranslationError,
)
from .models import TranslationResult

# Retry-eligible exception types (transient errors)
_RETRYABLE = (
    anthropic.APIConnectionError,
    anthropic.APITimeoutError,
    anthropic.RateLimitError,
    anthropic.InternalServerError,
)


async def _call_once(
    client: anthropic.AsyncAnthropic,
    messages: list[dict],
    label: str = "",
) -> anthropic.types.Message:
    """Single async Claude API call with jittered exponential backoff.

    Args:
        client: Async Anthropic client (shared across calls).
        messages: Conversation message history.
        label: Optional label for log messages (e.g., filename).
    """
    for attempt in range(MAX_RETRIES + 1):
        try:
            return await client.messages.create(
                model=MODEL,
                max_tokens=MAX_TOKENS,
                timeout=REQUEST_TIMEOUT,
                messages=messages,
                temperature=0,
            )
        except _RETRYABLE as e:
            if attempt == MAX_RETRIES:
                raise TranslationError(
                    f"[{label}] {type(e).__name__} after {MAX_RETRIES} retries: {e}"
                ) from e
            base_delay = BASE_DELAY * (2**attempt)
            jitter = random.uniform(0, base_delay * 0.5)
            delay = base_delay + jitter
            print(
                f"  [{label}] Retry {attempt + 1}/{MAX_RETRIES} "
                f"after {delay:.1f}s: {type(e).__name__}"
            )
            await asyncio.sleep(delay)
        except anthropic.APIError as e:
            raise TranslationError(f"[{label}] {type(e).__name__}: {e}") from e

    raise RuntimeError("Unreachable")  # pragma: no cover


async def call_claude_async(
    prompt: str,
    label: str = "",
    client: anthropic.AsyncAnthropic | None = None,
) -> TranslationResult:
    """Call Claude with automatic continuation for large responses.

    If the response hits max_tokens, the conversation is continued until
    the model emits an end_turn stop reason or MAX_CONTINUATIONS is reached.

    Args:
        prompt: The full prompt text.
        label: Identifier for log messages (e.g., relative file path).
        client: Optional shared async client; one is created if omitted.

    Returns:
        TranslationResult with the joined text and API usage metadata.
    """
    client = client or anthropic.AsyncAnthropic()
    messages: list[dict] = [{"role": "user", "content": prompt}]

    def _text(msg: anthropic.types.Message) -> str:
        if not msg.content or not hasattr(msg.content[0], "text"):
            raise TranslationError(
                f"[{label}] Unexpected API response: no text content block "
                f"(stop_reason={msg.stop_reason})"
            )
        return msg.content[0].text

    # Initial request
    message = await _call_once(client, messages, label)
    result_parts = [_text(message)]
    total_input = message.usage.input_tokens
    total_output = message.usage.output_tokens

    # Continuation loop for truncated responses
    continuations = 0
    while message.stop_reason == "max_tokens":
        continuations += 1
        if continuations > MAX_CONTINUATIONS:
            raise TranslationError(
                f"Response still incomplete after {MAX_CONTINUATIONS} continuations. "
                f"File may be too large to translate."
            )
        print(
            f"  [{label}] Response truncated, continuation "
            f"{continuations}/{MAX_CONTINUATIONS}..."
        )

        messages.append({"role": "assistant", "content": _text(message)})
        messages.append(
            {"role": "user", "content": "Continue exactly where you left off."}
        )

        message = await _call_once(client, messages, label)
        result_parts.append(_text(message))
        total_input += message.usage.input_tokens
        total_output += message.usage.output_tokens

    return TranslationResult(
        text="".join(result_parts).strip(),
        model=message.model,
        input_tokens=total_input,
        output_tokens=total_output,
        stop_reason=message.stop_reason,
        continuations=continuations,
    )
