<script>
	import { writable } from 'svelte/store';
	import { chatWithSpectrum } from '../services/api';

	export let hasSmilesPrediction = false;
	export let currentSmiles = null; // New prop to pass the current SMILES

	/**
	 * @typedef {Object} Message
	 * @property {number} id
	 * @property {string} avatar
	 * @property {string} name
	 * @property {string} message
	 * @property {string} timestamp
	 * @property {string} type
	 * @property {boolean} [thinking]
	 */

	/**
	 * @typedef {Object} ApiMessage
	 * @property {string} role
	 * @property {string} content
	 */

	/* message store ------------------------------------------------------- */
	/** @type {import('svelte/store').Writable<Message[]>} */
	const messages = writable([
		{
			id: 1,
			avatar: '/assets/images/spectra-avatar.png',
			name: 'Spectrum',
			message:
				"Hello! I'm Spectrum, your spectral-analysis assistant – ask me anything about SMILES or these spectra.",
			timestamp: new Date().toISOString(),
			type: 'assistant'
		}
	]);

	/* local state --------------------------------------------------------- */
	let userMessage = '';
	/** @type {HTMLDivElement} */
	let chatEl;
	let loading = false;
	let streamingMessageId = null;
	
	/* computed values ----------------------------------------------------- */
	$: canSend = userMessage.trim().length > 0 && !loading;

	/* helpers ------------------------------------------------------------- */
	const scrollToBottom = () =>
		setTimeout(() => chatEl?.scrollTo({ top: chatEl.scrollHeight, behavior: 'smooth' }), 40);

	$: scrollToBottom();               // whenever $messages changes
	
	// Handler for streaming message chunks
	function handleStreamChunk(chunk) {
		if (!streamingMessageId) return;
		
		messages.update(m => {
			const updatedMessages = [...m];
			const msgIndex = updatedMessages.findIndex(msg => msg.id === streamingMessageId);
			
			if (msgIndex !== -1) {
				updatedMessages[msgIndex] = {
					...updatedMessages[msgIndex],
					message: updatedMessages[msgIndex].message + chunk,
					thinking: false
				};
			}
			
			return updatedMessages;
		});
	}

	async function send() {
		if (!userMessage.trim() || loading) return;

		/* push user bubble */
		const userId = Date.now();
		messages.update((m) => [
			...m,
			{
				id: userId,
				avatar: '/assets/images/user-avatar.png',
				name: 'You',
				message: userMessage,
				timestamp: new Date().toISOString(),
				type: 'user'
			}
		]);

		/** @type {ApiMessage[]} */
		const history = $messages.map(({ type, message }) => ({ role: type, content: message }));
		userMessage = '';
		loading = true;

		/* Create assistant message for streaming */
		const assistantId = userId + 1;
		streamingMessageId = assistantId;
		
		messages.update((m) => [
			...m,
			{
				id: assistantId,
				avatar: '/assets/images/spectra-avatar.png',
				name: 'Spectrum',
				message: '', // Will be filled by streaming content
				timestamp: new Date().toISOString(),
				type: 'assistant',
				thinking: true
			}
		]);

		try {
			// Use streaming or regular response based on availability
			if ('ReadableStream' in window) {
				// Include current SMILES if available
				const chatOptions = {
					stream: true,
					onChunk: handleStreamChunk
				};
				
				if (currentSmiles) {
					chatOptions.smiles = currentSmiles;
				}
				
				await chatWithSpectrum(history, chatOptions);
				
				// Update message to show it's no longer streaming
				messages.update(m => {
					const updatedMessages = [...m];
					const msgIndex = updatedMessages.findIndex(msg => msg.id === streamingMessageId);
					
					if (msgIndex !== -1) {
						updatedMessages[msgIndex] = {
							...updatedMessages[msgIndex],
							thinking: false
						};
					}
					
					return updatedMessages;
				});
			} else {
				// Fall back to non-streaming if ReadableStream not supported
				const chatOptions = {};
				if (currentSmiles) {
					chatOptions.smiles = currentSmiles;
				}
				
				const { message } = await chatWithSpectrum(history, chatOptions);
				
				// Replace placeholder message with actual response
				messages.update(m => {
					const updatedMessages = [...m];
					const msgIndex = updatedMessages.findIndex(msg => msg.id === streamingMessageId);
					
					if (msgIndex !== -1) {
						updatedMessages[msgIndex] = {
							...updatedMessages[msgIndex],
							message,
							thinking: false
						};
					}
					
					return updatedMessages;
				});
			}
		} catch (e) {
			messages.update(m => {
				const updatedMessages = [...m];
				const msgIndex = updatedMessages.findIndex(msg => msg.id === streamingMessageId);
				
				if (msgIndex !== -1) {
					updatedMessages[msgIndex] = {
						...updatedMessages[msgIndex],
						message: 'Sorry – something went wrong. Try again?',
						thinking: false
					};
				}
				
				return updatedMessages;
			});
		} finally {
			loading = false;
			streamingMessageId = null;
		}
	}

	/**
	 * @param {KeyboardEvent} e
	 */
	function handleKey(e) {
		if (e.key === 'Enter' && !e.shiftKey) {
			e.preventDefault();
			send();
		}
	}
</script>

{#if hasSmilesPrediction}
	<div class="chat-window">
		<div class="messages" bind:this={chatEl}>
			{#each $messages as m}
				<div class="bubble {m.type} {m.thinking ? 'thinking' : ''}">
					<img class="avatar" alt={m.name} src={m.avatar} />
					<div class="content">
						<div class="header">
							<span class="name">{m.name}</span>
							<span class="time">{new Date(m.timestamp).toLocaleTimeString()}</span>
						</div>
						<p class="text">{m.message || 'Thinking...'}</p>
					</div>
				</div>
			{/each}
		</div>

		<div class="composer">
			<textarea
				rows="1"
				placeholder="Ask Spectrum…"
				bind:value={userMessage}
				on:keydown={handleKey}
				disabled={loading}
			/>

			<button
				class="send {canSend ? 'active' : ''} {loading ? 'loading' : ''}"
				disabled={!canSend}
				on:click={send}
				aria-label="Send message">
				<span class="send-icon">📤</span>
			</button>
		</div>
	</div>
{:else}
	<div class="placeholder">
		<span class="icon">💬</span>
		<p>Chat activates after your first prediction</p>
	</div>
{/if}

<style>
	:global(.chat-window) {
		display: flex;
		flex-direction: column;
		height: 100%;
	}

	.messages {
		flex: 1;
		overflow-y: auto;
		padding-right: .5rem;
		display: flex;
		flex-direction: column;
		gap: .75rem;
		margin-bottom: 1rem;
		scrollbar-width: thin;
	}

	.bubble {
		display: flex;
		gap: .6rem;
		max-width: 82%;
	}

	.bubble.user   { align-self: flex-end; flex-direction: row-reverse; }
	.bubble.assistant { align-self: flex-start; }

	.avatar {
		width: 26px; height: 26px; border-radius: 50%;
		flex-shrink: 0; background: var(--accent-soft);
	}

	.content {
		background: rgba(255,255,255,.75);
		border: 1px solid rgba(0,0,0,.05);
		border-radius: var(--enforce-pill);
		padding: .8rem 1.15rem;
		box-shadow: var(--shadow-sm);
		backdrop-filter: blur(8px);
	}

	.bubble.user .content {
		background: linear-gradient(135deg,var(--accent) 0%,var(--accent-secondary) 100%);
		color: #fff;
	}

	.header {
		display: flex; justify-content: space-between;
		margin-bottom: .35rem; font-size: .76rem;
		color: var(--text-secondary);
	}

	.bubble.user .header { color: rgba(255,255,255,.85); }

	.composer {
		display: flex; align-items: flex-end; gap: .65rem;
		border-top: 1px solid var(--surface-stroke);
		padding-top: .65rem;
		position: relative;
		z-index: 10;
	}

	textarea {
		flex: 1; resize: none;
		border: 1px solid var(--surface-stroke);
		border-radius: var(--enforce-pill);
		padding: .8rem 1rem;
		background: var(--surface-glass);
		font: 0.9rem/1.4 system-ui;
		position: relative;
		z-index: 15;
	}

	textarea:focus-visible { outline: 2px solid var(--accent); }

	/* default (disabled) – soft grey pill */
	button.send {
		width: 46px; height: 46px;
		border: none;
		border-radius: 50%;
		background: var(--accent-soft);
		color: var(--text-tertiary);
		cursor: not-allowed;
		box-shadow: none;
		transition: transform .15s ease, background .15s ease;
		font-size: 1.25rem;
		display: flex;
		align-items: center;
		justify-content: center;
		position: relative;
		z-index: 20;
	}

	/* when the user can send – bright purple gradient */
	button.send.active {
		cursor: pointer;
		background: linear-gradient(135deg, var(--accent) 0%, var(--accent-secondary) 100%);
		color: #fff;
		box-shadow: 0 4px 12px rgba(120,121,255,.35);
	}
	button.send.active:hover { transform: translateY(-2px); }

	/* mini-loader while waiting for Spectrum's reply */
	button.send.loading {
		background: var(--accent-soft);
		cursor: progress;
	}

	.send-icon { 
		position: relative;
		z-index: 2;      /* sit on top of the gradient */
		line-height: 1;   /* no extra vertical space */
	}

	/* Add animation for "thinking" state */
	.bubble.thinking .text {
		opacity: 0.7;
	}
	
	.bubble.thinking .content {
		animation: pulse 1.5s infinite;
	}
	
	@keyframes pulse {
		0% { opacity: 0.7; }
		50% { opacity: 1; }
		100% { opacity: 0.7; }
	}
</style> 