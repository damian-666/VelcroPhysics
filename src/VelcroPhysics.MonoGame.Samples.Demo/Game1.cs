﻿using System;
using System.Collections.Generic;
using System.Reflection;
using Genbox.VelcroPhysics.MonoGame.Samples.Demo.MediaSystem;
using Genbox.VelcroPhysics.MonoGame.Samples.Demo.MediaSystem.Graphics;
using Genbox.VelcroPhysics.MonoGame.Samples.Demo.ScreenSystem;
using Genbox.VelcroPhysics.Utilities;
using Microsoft.Xna.Framework;
using Microsoft.Xna.Framework.Graphics;
using Microsoft.Xna.Framework.Input;

namespace Genbox.VelcroPhysics.MonoGame.Samples.Demo
{
    /// <summary>
    /// This is the main type for the samples
    /// </summary>
    public sealed class Game1 : Game
    {
        private readonly GraphicsDeviceManager _graphics;
        private readonly List<GameScreen> _screens = new List<GameScreen>();
        private readonly List<GameScreen> _screensToUpdate = new List<GameScreen>();
        private readonly List<RenderTarget2D> _transitions = new List<RenderTarget2D>();
        private InputHelper _input;

        private bool _isExiting;
        private LineBatch _lineBatch;
        private MenuScreen _menuScreen;
        private QuadRenderer _quadRenderer;
        private SpriteBatch _spriteBatch;

#if WINDOWS
        private FrameRateCounter _counter;
        private bool _showFps;
#endif

        public Game1()
        {
            Window.Title = "Velcro Physics Samples";
            _graphics = new GraphicsDeviceManager(this);
            ConvertUnits.SetDisplayUnitToSimUnitRatio(24f);
            IsFixedTimeStep = true;
            Content.RootDirectory = "Content";
        }

        protected override void Initialize()
        {
            _graphics.PreferMultiSampling = true;
            _graphics.PreferredBackBufferWidth = GraphicsDevice.DisplayMode.Width - 80;
            _graphics.PreferredBackBufferHeight = GraphicsDevice.DisplayMode.Height - 120;

#if WINDOWS
            _graphics.IsFullScreen = false;
#elif XBOX
            _graphics.IsFullScreen = true;
#endif

            _graphics.ApplyChanges();

            _input = new InputHelper();
            _isExiting = false;

#if WINDOWS
            _counter = new FrameRateCounter();
            _showFps = false;
#endif

            base.Initialize();
        }

        protected override void LoadContent()
        {
            base.LoadContent();

            Managers.TextureManager = new TextureManager(Content, GraphicsDevice);
            Managers.FontManager = new FontManager(Content);
            Managers.SoundManager = new SoundManager(Content);

            _spriteBatch = new SpriteBatch(GraphicsDevice);
            _lineBatch = new LineBatch(GraphicsDevice);
            _quadRenderer = new QuadRenderer(GraphicsDevice);

            _input.LoadContent(GraphicsDevice.Viewport);

#if WINDOWS
            _counter.LoadContent();
#endif

            // Create rendertarget for transitions
            PresentationParameters pp = GraphicsDevice.PresentationParameters;
            _transitions.Add(new RenderTarget2D(GraphicsDevice, pp.BackBufferWidth, pp.BackBufferHeight, false, SurfaceFormat.Color, pp.DepthStencilFormat, pp.MultiSampleCount, RenderTargetUsage.DiscardContents));

            _menuScreen = new MenuScreen();

            List<Type> demosToLoad = new List<Type>();
            Assembly samplesFramework = Assembly.GetExecutingAssembly();
            foreach (Type sampleType in samplesFramework.GetTypes())
            {
                if (sampleType.IsSubclassOf(typeof(PhysicsDemoScreen)))
                    demosToLoad.Add(sampleType);
            }
            demosToLoad.Add(demosToLoad[0]); // HACK: Load the first sample two times, since some delayed creation stuff with the rendertargets always breaks the first preview picture.

            bool firstPreview = true;
            foreach (Type sampleType in demosToLoad)
            {
                PhysicsDemoScreen demoScreen = (PhysicsDemoScreen)samplesFramework.CreateInstance(sampleType.ToString());

#if WINDOWS
                if (!firstPreview)
                    Console.WriteLine("Loading demo: " + demoScreen.GetTitle());
#endif

                demoScreen.Framework = this;
                demoScreen.IsExiting = false;

                demoScreen.Sprites = _spriteBatch;
                demoScreen.Lines = _lineBatch;
                demoScreen.Quads = _quadRenderer;

                demoScreen.LoadContent();

                // "Abuse" transition rendertarget to render screen preview
                GraphicsDevice.SetRenderTarget(_transitions[0]);
                GraphicsDevice.Clear(Color.Transparent);

                _quadRenderer.Begin();
                _quadRenderer.Render(Vector2.Zero, new Vector2(_transitions[0].Width, _transitions[0].Height), null, true, Colors.Grey, Color.White * 0.3f);
                _quadRenderer.End();

                // Update ensures that the screen is fully visible, we "cover" it so that no physics are run
                demoScreen.Update(new GameTime(demoScreen.TransitionOnTime, demoScreen.TransitionOnTime), true, false);
                demoScreen.Draw(new GameTime());
                demoScreen.Draw(new GameTime());

                using (RenderTarget2D preview = new RenderTarget2D(GraphicsDevice, pp.BackBufferWidth / 2, pp.BackBufferHeight / 2, false, SurfaceFormat.Color, pp.DepthStencilFormat, pp.MultiSampleCount, RenderTargetUsage.DiscardContents))
                {
                    GraphicsDevice.SetRenderTarget(preview);
                    GraphicsDevice.Clear(Color.Transparent);

                    _spriteBatch.Begin();
                    _spriteBatch.Draw(_transitions[0], preview.Bounds, Color.White);
                    _spriteBatch.End();

                    GraphicsDevice.SetRenderTarget(null);

                    demoScreen.ExitScreen();
                    demoScreen.Update(new GameTime(demoScreen.TransitionOffTime, demoScreen.TransitionOffTime), true, false);

                    if (!firstPreview)
                        _menuScreen.AddMenuItem(demoScreen, preview);
                    else
                        firstPreview = false;
                }
            }

            AddScreen(new BackgroundScreen());
            AddScreen(_menuScreen);

            //TODO: Can't call this in MonoGame at the moment
            //ResetElapsedTime();
        }

        protected override void UnloadContent()
        {
            foreach (GameScreen screen in _screens)
            {
                screen.UnloadContent();
            }
            base.UnloadContent();
        }

        /// <summary>
        /// Allows each screen to run logic.
        /// </summary>
        protected override void Update(GameTime gameTime)
        {
            // Read the keyboard and gamepad.
            _input.Update(gameTime);

#if WINDOWS
            // Update the framerate counter
            _counter.Update(gameTime);
#endif

            if ((_input.IsNewButtonPress(Buttons.Y) || _input.IsNewKeyPress(Keys.F5)) && !(_screens[_screens.Count - 1] is OptionsScreen))
                AddScreen(new OptionsScreen());

#if WINDOWS
            if (_input.IsNewKeyPress(Keys.F11))
                _showFps = !_showFps;

            if (_input.IsNewKeyPress(Keys.F12))
                _graphics.ToggleFullScreen();
#endif

            // Make a copy of the master screen list, to avoid confusion if
            // the process of updating one screen adds or removes others.
            _screensToUpdate.Clear();
            _screensToUpdate.AddRange(_screens);

            bool otherScreenHasFocus = !IsActive;
            bool coveredByOtherScreen = false;

            // Loop as long as there are screens waiting to be updated.
            while (_screensToUpdate.Count > 0)
            {
                // Pop the topmost screen off the waiting list.
                GameScreen screen = _screensToUpdate[_screensToUpdate.Count - 1];

                _screensToUpdate.RemoveAt(_screensToUpdate.Count - 1);

                // Update the screen.
                screen.Update(gameTime, otherScreenHasFocus, coveredByOtherScreen);

                if (screen.ScreenState == ScreenState.TransitionOn || screen.ScreenState == ScreenState.Active)
                {
                    // If this is the first active screen we came across,
                    // give it a chance to handle input.
                    if (!otherScreenHasFocus && !_isExiting)
                    {
                        _input.ShowCursor = screen.HasCursor;
                        screen.HandleInput(_input, gameTime);
                        otherScreenHasFocus = true;
                    }

                    // If this is an active non-popup, inform any subsequent
                    // screens that they are covered by it.
                    if (!screen.IsPopup)
                        coveredByOtherScreen = true;
                }
            }

            if (_isExiting && _screens.Count == 0)
                Exit();

            base.Update(gameTime);
        }

        protected override void Draw(GameTime gameTime)
        {
            int transitionCount = 0;
            foreach (GameScreen screen in _screens)
            {
                if (screen.ScreenState == ScreenState.TransitionOn || screen.ScreenState == ScreenState.TransitionOff)
                {
                    transitionCount++;
                    if (_transitions.Count < transitionCount)
                    {
                        PresentationParameters pp = GraphicsDevice.PresentationParameters;
                        _transitions.Add(new RenderTarget2D(GraphicsDevice, pp.BackBufferWidth, pp.BackBufferHeight, false, SurfaceFormat.Color, pp.DepthStencilFormat, pp.MultiSampleCount, RenderTargetUsage.DiscardContents));
                    }
                    GraphicsDevice.SetRenderTarget(_transitions[transitionCount - 1]);
                    GraphicsDevice.Clear(Color.Transparent);
                    screen.Draw(gameTime);
                    GraphicsDevice.SetRenderTarget(null);
                }
            }

            GraphicsDevice.Clear(Color.Black);

            transitionCount = 0;
            foreach (GameScreen screen in _screens)
            {
                if (screen.ScreenState == ScreenState.Hidden)
                    continue;

                if (screen.ScreenState == ScreenState.TransitionOn || screen.ScreenState == ScreenState.TransitionOff)
                {
                    _spriteBatch.Begin(0, BlendState.AlphaBlend);
                    if (screen is PhysicsDemoScreen)
                    {
                        Vector2 position = Vector2.Lerp(_menuScreen.PreviewPosition, new Vector2(GraphicsDevice.Viewport.Width, GraphicsDevice.Viewport.Height) / 2f, 1f - screen.TransitionPosition);
                        _spriteBatch.Draw(_transitions[transitionCount], position, null, Color.White * Math.Min(screen.TransitionAlpha / 0.2f, 1f), 0f,
                            new Vector2(_transitions[transitionCount].Width, _transitions[transitionCount].Height) / 2f, 0.5f + 0.5f * (1f - screen.TransitionPosition), SpriteEffects.None, 0f);
                    }
                    else
                        _spriteBatch.Draw(_transitions[transitionCount], Vector2.Zero, Color.White * screen.TransitionAlpha);

                    _spriteBatch.End();

                    transitionCount++;
                }
                else
                    screen.Draw(gameTime);
            }

            _input.Draw(_spriteBatch);

#if WINDOWS
            if (_showFps)
                _counter.Draw(_spriteBatch);
#endif
            base.Draw(gameTime);
        }

        public void ExitGame()
        {
            foreach (GameScreen screen in _screens)
            {
                screen.ExitScreen();
            }
            _isExiting = true;
        }

        /// <summary>
        /// Adds a new screen to the screen manager.
        /// </summary>
        public void AddScreen(GameScreen screen)
        {
            screen.Framework = this;
            screen.IsExiting = false;

            screen.Sprites = _spriteBatch;
            screen.Lines = _lineBatch;
            screen.Quads = _quadRenderer;

            // Tell the screen to load content.
            screen.LoadContent();

            // Loading my take a while so elapsed time is reset to prevent hick-ups
            //TODO: Can't call this in MonoGame at the moment
            //ResetElapsedTime();

            _screens.Add(screen);
        }

        /// <summary>
        /// Removes a screen from the screen manager. You should normally
        /// use GameScreen.ExitScreen instead of calling this directly, so
        /// the screen can gradually transition off rather than just being
        /// instantly removed.
        /// </summary>
        public void RemoveScreen(GameScreen screen)
        {
            // Tell the screen to unload content.
            screen.UnloadContent();
            _screens.Remove(screen);
            _screensToUpdate.Remove(screen);
        }

        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                _graphics.Dispose();
                _lineBatch.Dispose();
                _quadRenderer.Dispose();
                _spriteBatch.Dispose();
            }

            base.Dispose(disposing);
        }
    }
}